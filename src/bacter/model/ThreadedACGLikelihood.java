package bacter.model;

import bacter.*;
import beast.app.BeastMCMC;
import beast.core.Description;
import beast.core.Input;
import beast.core.State;
import beast.core.util.Log;
import beast.evolution.alignment.Alignment;
import beast.evolution.branchratemodel.BranchRateModel;
import beast.evolution.branchratemodel.StrictClockModel;
import beast.evolution.likelihood.*;
import beast.evolution.sitemodel.SiteModel;
import beast.evolution.substitutionmodel.SubstitutionModel;
import beast.evolution.tree.Node;
import com.google.common.collect.LinkedHashMultiset;
import com.google.common.collect.Multiset;

import java.util.*;
import java.util.concurrent.*;


/**
 * @author Tim Vaughan <tgvaughan@gmail.com>
 */
@Description("Probability of sequence data given recombination graph.")
public class ThreadedACGLikelihood extends GenericTreeLikelihood {

    final public Input<Integer> maxNrOfThreadsInput = new Input<>("threads","maximum number of threads to use, if less than 1 the number of threads in BeastMCMC is used (default -1)", -1);

    public Input<Locus> th_locusInput = new Input<>(
            "locus",
            "Locus associated with alignment to evaluate probability of.",
            Input.Validate.REQUIRED);

    public Input<Boolean> th_useAmbiguitiesInput = new Input<>(
            "useAmbiguities",
            "Whether sites containing ambiguous states should be handled " +
                    "instead of ignored (the default)", false);

    protected ConversionGraph acg;
    protected SiteModel.Base siteModel;
    protected BranchRateModel.Base branchRateModel;
    protected SubstitutionModel.Base substitutionModel;
    protected Alignment alignment;
    protected Locus locus;
    protected int nStates;

    protected ConcurrentMap<Region, Multiset<int[]>> patterns;
    protected ConcurrentMap<Region, Multiset<int[]>> storedPatterns;
    protected ConcurrentMap<Region, double[]> patternLogLikelihoods;
    protected ConcurrentMap<Region, double[]> storedPatternLogLikelihoods;
    protected ConcurrentMap<Region, double[]> rootPartials;
    protected ConcurrentMap<Region, double[]> storedRootPartials;
    protected ConcurrentMap<Region, List<Integer>> constantPatterns;
    protected ConcurrentMap<Region, List<Integer>> storedConstantPatterns;
    protected ConcurrentMap<Region, LikelihoodCore> likelihoodCores;
    protected ConcurrentMap<Region, LikelihoodCore> storedLikelihoodCores;
    protected ConcurrentMap<Region, Double> regionLogLikelihoods;
    protected ConcurrentMap<Region, Double> storedRegionLogLikelihoods;

    /**
     * Memory for transition probabilities.
     */
    protected double [][] probabilities;

    public ThreadedACGLikelihood() {
        // We allow alignments to be specified using Locus objects.
        dataInput.setRule(Input.Validate.OPTIONAL);
    }


    /** calculation engine **/
    private ExecutorService pool = null;
    private final List<Callable<Double>> likelihoodCallers = new ArrayList<>();


    /** number of threads to use, changes when threading causes problems **/
    private int threadCount;
    public double[] logPByThread;



    /**
     * Cached transition probabilities for CF edges.
     */
    double [][][][] cfTransitionProbs;
    int cacheHits = 0;
    int [] cacheHitsByThread;
    int cacheMisses = 0;
    int [] cacheMissesByThread;

    MarginalTree[] marginalTreePerRegion;
    List<Region> currentRegionList;

    @Override
    public void initAndValidate() {

        if (treeInput.get() instanceof ConversionGraph)
            acg = (ConversionGraph)treeInput.get();
        else
            throw new IllegalArgumentException("'Tree' input to ACGLikelihood must " +
                    "be of type ConversionGraph.");

        threadCount = BeastMCMC.m_nThreads;
        if (maxNrOfThreadsInput.get() > 0) {
            threadCount = Math.min(maxNrOfThreadsInput.get(), BeastMCMC.m_nThreads);
        }
        String instanceCount = System.getProperty("beast.instance.count");
        if (instanceCount != null && instanceCount.length() > 0) {
            threadCount = Integer.parseInt(instanceCount);
        }

        pool = Executors.newFixedThreadPool(threadCount);
        logPByThread = new double[threadCount];

        locus = th_locusInput.get();
        if (locus.hasAlignment()) {
            alignment = locus.getAlignment();
        } else {
            if (dataInput.get() != null)
                alignment = dataInput.get();
            else
                throw new IllegalArgumentException("No alignment associated with " +
                        "locus " + locus.getID() + " provided to ACGLikelihood " +
                        "and none given explicitly.");
        }

        nStates = alignment.getMaxStateCount();

        siteModel = (SiteModel.Base) siteModelInput.get();
        substitutionModel = (SubstitutionModel.Base) siteModel.getSubstitutionModel();

        if (branchRateModelInput.get() != null) {
            branchRateModel = branchRateModelInput.get();

            if (!(branchRateModel instanceof StrictClockModel))
                throw new IllegalArgumentException("ACGLikelihood currently only" +
                        "supports strict clock models.");
        } else
            branchRateModel = new StrictClockModel();

        patterns = new ConcurrentHashMap<>();
        storedPatterns = new ConcurrentHashMap<>();
        patternLogLikelihoods = new ConcurrentHashMap<>();
        storedPatternLogLikelihoods = new ConcurrentHashMap<>();
        rootPartials = new ConcurrentHashMap<>();
        storedRootPartials = new ConcurrentHashMap<>();
        constantPatterns = new ConcurrentHashMap<>();
        storedConstantPatterns = new ConcurrentHashMap<>();
        likelihoodCores = new ConcurrentHashMap<>();
        storedLikelihoodCores = new ConcurrentHashMap<>();
        regionLogLikelihoods = new ConcurrentHashMap<>();
        storedRegionLogLikelihoods = new ConcurrentHashMap<>();


        // Allocate transition probability memory:
        // (Only the first nStates*nStates elements are usually used.)
        probabilities = new double[threadCount][(nStates+1)*(nStates+1)];

        stack = new ArrayDeque[threadCount];
        for(int i=0; i<threadCount;i++) stack[i] = new ArrayDeque<>();
        postOrderNodes = new MarginalNode[threadCount][];
        cacheHitsByThread = new int[threadCount];
        cacheMissesByThread = new int[threadCount];

    }

    protected double scaleFactor = 1.0;
    protected double scaleFactorMultiplier = 1.01;
    protected double maxScaleFactor = 100.0;
    int numberOfUsedThreads;

    @Override
    public double calculateLogP() {

        numberOfUsedThreads = Math.min(threadCount,acg.getRegions(locus).size());

        while (true) {
            doLogPCalculation();

            if (logP>Double.NEGATIVE_INFINITY || scaleFactor>=maxScaleFactor) {
                return logP;
            }

            scaleFactor *= scaleFactorMultiplier;
            Log.warning.println("Turning on scaling to prevent numeric instability. Scale factor: " + scaleFactor);

            for (LikelihoodCore core : likelihoodCores.values())
                core.setUseScaling(scaleFactor);

            requiresRecalculation();
        }
    }



    protected void doLogPCalculation() {
        updatePatterns();
        updateCores();
        preComputeCFTransitionProbs();

        for(int i = 0; i< threadCount; i++) {
            logPByThread[i] = 0;
            cacheHitsByThread[i] = 0;
            cacheMissesByThread[i] = 0;
        }

        logP = 0.0;
        regionLogLikelihoods.keySet().retainAll(acg.getRegions(locus));

        if(numberOfUsedThreads <= 1) {

            for (Region region : acg.getRegions(locus)) {

                if (!regionLogLikelihoods.containsKey(region)) {
                    traverseNoRecurse(new MarginalTree(acg, region.activeConversions).getRoot(), region, 0);

                    double regionLogP = 0.0;
                    int i = 0;
                    for (int[] pattern : patterns.get(region).elementSet()) {
                        regionLogP += patternLogLikelihoods.get(region)[i]
                                * patterns.get(region).count(pattern);
                        i += 1;
                    }
                    regionLogLikelihoods.put(region, regionLogP);

                    logP += regionLogP;
                } else {
                    logP += regionLogLikelihoods.get(region);
                }
            }
        }
        else{
            prepareLogPCalculationInThreads();
        }
    }

    protected void prepareLogPCalculationInThreads(){

        int regionCount = acg.getRegions(locus).size();

        /* regionsPerThread is the range of regions for each thread, except the last one*/
        int regionsPerThread = Math.round(regionCount / threadCount);
        int regionNumber = 0;


        likelihoodCallers.clear();
        currentRegionList = new ArrayList<>();
        marginalTreePerRegion = new MarginalTree[acg.getRegions(locus).size()];
        currentRegionList.addAll(acg.getRegions(locus));

        /*Setup the region list (for first threads) and create caller. The last thread is set up separately*/
        for (int i = 0; i < numberOfUsedThreads - 1; i++) {
            List<Region> regionListForThread = new ArrayList<>();

            for (int j = 0; j < regionsPerThread; j++) {
                regionListForThread.add(currentRegionList.get(regionNumber));
                marginalTreePerRegion[regionNumber] = new MarginalTree(acg, acg.getRegions(locus).get(regionNumber).activeConversions);
                regionNumber += 1;
            }
            if(regionListForThread.size()>0) {
                likelihoodCallers.add(new ThreadedACGLikelihood.ACGLikelihoodCaller(regionListForThread, i));
            }
        }

        /*last thread might have different number of regions*/
        List<Region> regionListForThread = new ArrayList<>();
        while (regionNumber < regionCount) {
            regionListForThread.add(currentRegionList.get(regionNumber));
            marginalTreePerRegion[regionNumber] = new MarginalTree(acg, acg.getRegions(locus).get(regionNumber).activeConversions);
            regionNumber += 1;
        }
        if(regionListForThread.size()>0) {
            likelihoodCallers.add(new ThreadedACGLikelihood.ACGLikelihoodCaller(regionListForThread, numberOfUsedThreads -1));
        }

        /*perform calculation in threads and sum over the logP Array*/
        calculateThreadedLogP();

    }

    private void calculateThreadedLogP() {
        try {

            pool.invokeAll(likelihoodCallers);

            logP = 0.0;
            for(int i = 0; i< numberOfUsedThreads; i++){

                logP += logPByThread[i];
                cacheHits += cacheHitsByThread[i];
                cacheMisses += cacheMissesByThread[i];

            }


        } catch (RejectedExecutionException | InterruptedException e) {
            e.printStackTrace();
            System.exit(0);
        }
    }

    class ACGLikelihoodCaller implements Callable<Double> {
        private final List<Region> regionListForThread;
        private final int threadNr;

        public ACGLikelihoodCaller(List<Region> regionListForThread, int threadNr) {
            this.regionListForThread = regionListForThread;
            this.threadNr = threadNr;
        }

        public Double call() throws Exception {
            try {

                for (Region region : regionListForThread) {

                    if (!regionLogLikelihoods.containsKey(region)) {
                        traverseNoRecurse(marginalTreePerRegion[currentRegionList.indexOf(region)].getRoot(), region, threadNr);

                        double regionLogP = 0.0;
                        int i = 0;
                        for (int[] pattern : patterns.get(region).elementSet()) {
                            regionLogP += patternLogLikelihoods.get(region)[i]
                                    * patterns.get(region).count(pattern);
                            i += 1;
                        }
                        regionLogLikelihoods.put(region, regionLogP);
                        logPByThread[threadNr] += regionLogP;

                    } else {
                        logPByThread[threadNr] += regionLogLikelihoods.get(region);
                    }

                }

            } catch (Exception e) {
                System.err.println("Something went wrong in thread " + threadNr);
                e.printStackTrace();
                System.exit(0);
            }
            return logPByThread[threadNr];
        }

    }

    /**
     * Ensure pattern counts are up to date.
     */
    protected void updatePatterns() {
        List<Region> regionList = acg.getRegions(locus);

        // Remove stale pattern sets
        patterns.keySet().retainAll(regionList);
        patternLogLikelihoods.keySet().retainAll(regionList);
        rootPartials.keySet().retainAll(regionList);
        constantPatterns.keySet().retainAll(regionList);

        for (Region region : regionList) {

            if (patterns.containsKey(region))
                continue;

            // Add new pattern set
            Multiset<int[]> patSet = LinkedHashMultiset.create();
            for (int j=region.leftBoundary; j<region.rightBoundary; j++) {
                int [] pat = alignment.getPattern(alignment.getPatternIndex(j));
                patSet.add(pat);
            }

            //todo: check adjustment (circular genome)
            if (region.leftBoundary > region.rightBoundary) {
                for (int j=region.leftBoundary; j<acg.getTotalConvertibleSequenceLength(); j++) {
                    int[] pat = alignment.getPattern(alignment.getPatternIndex(j));
                    patSet.add(pat);
                }
                for (int j=0; j<region.rightBoundary; j++) {
                    int[] pat = alignment.getPattern(alignment.getPatternIndex(j));
                    patSet.add(pat);
                }
            }

            patterns.put(region, patSet);

            // Allocate memory for corresponding log likelihoods and root partials
            patternLogLikelihoods.put(region, new double[patSet.elementSet().size()]);
            rootPartials.put(region, new double[patSet.elementSet().size()*nStates]);

            // Compute corresponding constant pattern list
            List<Integer> constantPatternList = new ArrayList<>();

            int patternIdx = 0;
            for (int[] pattern : patSet.elementSet()) {
                boolean isConstant = true;
                for (int i=1; i<pattern.length; i++)
                    if (pattern[i] != pattern[0]) {
                        isConstant = false;
                        break;
                    }

                if (isConstant) {
                    if (alignment.getDataType().isAmbiguousCode(pattern[0])) {
                        if (th_useAmbiguitiesInput.get()) {
                            for (int state : alignment.getDataType().getStatesForCode(pattern[0]))
                                constantPatternList.add(patternIdx * nStates + state);
                        }
                    } else {
                        constantPatternList.add(patternIdx * nStates + pattern[0]);
                    }
                }

                patternIdx += 1;
            }

            constantPatterns.put(region, constantPatternList);
        }
    }


    /**
     * Initialize likelihood cores.
     */
    protected void updateCores() {

        List<Region> regionList = acg.getRegions(locus);
        likelihoodCores.keySet().retainAll(regionList);

        for (Region region : regionList) {

            if (likelihoodCores.containsKey(region))
                continue;

            LikelihoodCore likelihoodCore;
            if (nStates==4)
                likelihoodCore = new BeerLikelihoodCore4();
            else
                likelihoodCore = new BeerLikelihoodCore(nStates);

            likelihoodCores.put(region, likelihoodCore);

            likelihoodCore.initialize(acg.getNodeCount(),
                    patterns.get(region).elementSet().size(),
                    siteModel.getCategoryCount(),
                    true, th_useAmbiguitiesInput.get());

            if (scaleFactor>1.0)
                likelihoodCore.setUseScaling(scaleFactor);

            if (th_useAmbiguitiesInput.get())
                setPartials(likelihoodCore, patterns.get(region));
            else
                setStates(likelihoodCore, patterns.get(region));

            int intNodeCount = acg.getNodeCount()/2;
            for (int i=0; i<intNodeCount; i++)
                likelihoodCore.createNodePartials(intNodeCount+1+i);
        }
    }


    /**
     * Set leaf states in a likelihood core.
     *
     * @param lhc       likelihood core object
     * @param patterns  leaf state patterns
     */
    void setStates(LikelihoodCore lhc, Multiset<int[]> patterns) {

        for (Node node : acg.getExternalNodes()) {
            int[] states = new int[patterns.elementSet().size()];
            int taxon = alignment.getTaxonIndex(node.getID());
            int i=0;
            for (int [] pattern : patterns.elementSet()) {
                int code = pattern[taxon];
                int[] statesForCode = alignment.getDataType().getStatesForCode(code);
                if (statesForCode.length==1)
                    states[i] = statesForCode[0];
                else
                    states[i] = code; // Causes ambiguous states to be ignored.

                i += 1;
            }
            lhc.setNodeStates(node.getNr(), states);
        }
    }


    /**
     * Set leaf partials in likelihood core.
     *
     * @param lhc likelihood core object
     * @param patterns leaf state patterns
     */
    protected void setPartials(LikelihoodCore lhc, Multiset<int[]> patterns) {
        for (Node node : acg.getExternalNodes()) {
            int nStates = alignment.getDataType().getStateCount();
            double[] partials = new double[patterns.elementSet().size() * nStates];
            int k = 0;
            int iTaxon = alignment.getTaxonIndex(node.getID());
            for (int[] pattern : patterns.elementSet()) {
                int code = pattern[iTaxon];
                boolean[] stateSet = alignment.getDataType().getStateSet(code);
                for (int iState = 0; iState < nStates; iState++) {
                    partials[k++] = (stateSet[iState] ? 1.0 : 0.0);
                }
            }
            lhc.setNodePartials(node.getNr(), partials);
        }
    }


    /**
     * Pre-compute transition probabilities for CF edges.
     */
    void preComputeCFTransitionProbs() {
        if (cfTransitionProbs == null)
            cfTransitionProbs = new double[threadCount][acg.getNodeCount()-1][siteModel.getCategoryCount()][(nStates+1)*(nStates+1)];

        for(int th = 0; th< threadCount; th++) {
            for (int ni = 0; ni < acg.getNodeCount() - 1; ni++) {
                Node node = acg.getNode(ni);
                for (int ci = 0; ci < siteModel.getCategoryCount(); ci++) {
                    double jointBranchRate = siteModel.getRateForCategory(ci, node)
                            * branchRateModel.getRateForBranch(node);
                    double parentHeight = node.getParent().getHeight();
                    double nodeHeight = node.getHeight();

                    substitutionModel.getTransitionProbabilities(
                            node,
                            parentHeight,
                            nodeHeight,
                            jointBranchRate,
                            cfTransitionProbs[th][ni][ci]);
                }
            }
        }
    }

    ArrayDeque<MarginalNode>[] stack;
    MarginalNode[][] postOrderNodes;


    /**
     * Assemble list of marginal tree nodes for post-order traversal.
     *
     * @param root root of marginal tree
     */
    @SuppressWarnings("deprecation")
    void computePostOrder(MarginalNode root, int threadNumber) {

        if (postOrderNodes[threadNumber] == null)
            postOrderNodes[threadNumber] = new MarginalNode[acg.getNodeCount()];

        stack[threadNumber].clear();
        int i = 0;

        stack[threadNumber].push(root);
        while (!stack[threadNumber].isEmpty()) {
            MarginalNode n = stack[threadNumber].pop();
            if (!n.isLeaf()) {
                stack[threadNumber].push((MarginalNode)n.getLeft());
                stack[threadNumber].push((MarginalNode)n.getRight());
            }
            postOrderNodes[threadNumber][postOrderNodes[threadNumber].length-1-(i++)] = n;
        }
    }

    /**
     * Traverse a marginal tree, computing partial likelihoods on the way.
     * This version avoids potentially-expensive recursive function calls.
     *
     * @param root Tree node
     * @param region region
     */
    void traverseNoRecurse(MarginalNode root, Region region, int threadNumber) {

        computePostOrder(root, threadNumber);

        LikelihoodCore lhc = likelihoodCores.get(region);

        for (MarginalNode node : postOrderNodes[threadNumber]) {

            if (!node.isRoot()) {
                lhc.setNodeMatrixForUpdate(node.getNr());

                boolean cfEdge = node.cfNodeNr>=0
                        && !acg.getNode(node.cfNodeNr).isRoot()
                        && acg.getNode(node.cfNodeNr).getParent().getNr()
                        == ((MarginalNode)node.getParent()).cfNodeNr;

                if (!cfEdge) {
                    cacheMissesByThread[threadNumber] += 1;

                    for (int i = 0; i < siteModel.getCategoryCount(); i++) {
                        double jointBranchRate = siteModel.getRateForCategory(i, node)
                                * branchRateModel.getRateForBranch(node);
                        double parentHeight = node.getParent().getHeight();
                        double nodeHeight = node.getHeight();

                        substitutionModel.getTransitionProbabilities(
                                node,
                                parentHeight,
                                nodeHeight,
                                jointBranchRate,
                                probabilities[threadNumber]);
                        lhc.setNodeMatrix(node.getNr(), i, probabilities[threadNumber]);
                    }
                } else {
                    cacheHitsByThread[threadNumber] += 1;

                    for (int i=0; i<siteModel.getCategoryCount(); i++) {
                        lhc.setNodeMatrix(node.getNr(), i, cfTransitionProbs[threadNumber][node.cfNodeNr][i]);
                    }
                }
            }

            if (!node.isLeaf()) {

                // LikelihoodCore only supports binary trees.
                List<Node> children = node.getChildren();
                lhc.setNodePartialsForUpdate(node.getNr());
                lhc.setNodeStatesForUpdate(node.getNr());
                lhc.calculatePartials(children.get(0).getNr(),
                        children.get(1).getNr(), node.getNr());

                if (node.isRoot()) {
                    double[] frequencies = substitutionModel.getFrequencies();
                    double[] proportions = siteModel.getCategoryProportions(node);
                    lhc.integratePartials(node.getNr(), proportions,
                            rootPartials.get(region));

                    for (int idx : constantPatterns.get(region)) {
                        rootPartials.get(region)[idx]
                                += siteModel.getProportionInvariant();
                    }

                    lhc.calculateLogLikelihoods(rootPartials.get(region),
                            frequencies, patternLogLikelihoods.get(region));
                }
            }

        }
    }

    @Override
    public List<String> getArguments() {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    @Override
    public List<String> getConditions() {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    @Override
    public void sample(State state, Random random) {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    @Override
    protected boolean requiresRecalculation() {

        if (acg.clonalFrameIsDirty()
                || siteModel.isDirtyCalculation()
                || branchRateModel.isDirtyCalculation())
            regionLogLikelihoods.clear();

        return true;
    }

    @Override
    public void store() {
        storedPatterns.clear();
        storedPatterns.putAll(patterns);

        storedPatternLogLikelihoods.clear();
        storedPatternLogLikelihoods.putAll(patternLogLikelihoods);

        storedConstantPatterns.clear();
        storedConstantPatterns.putAll(constantPatterns);

        storedRootPartials.clear();
        storedRootPartials.putAll(rootPartials);

        storedLikelihoodCores.clear();
        storedLikelihoodCores.putAll(likelihoodCores);

        storedRegionLogLikelihoods.clear();
        storedRegionLogLikelihoods.putAll(regionLogLikelihoods);

        super.store();
    }

    @Override
    public void restore() {
        ConcurrentMap<Region, Multiset<int[]>> tmpPatterns = patterns;
        patterns = storedPatterns;
        storedPatterns = tmpPatterns;

        ConcurrentMap<Region, double[]> tmpPatternLogLikelihoods = patternLogLikelihoods;
        patternLogLikelihoods = storedPatternLogLikelihoods;
        storedPatternLogLikelihoods = tmpPatternLogLikelihoods;

        ConcurrentMap<Region, double[]> tmpRootPartials = rootPartials;
        rootPartials = storedRootPartials;
        storedRootPartials = tmpRootPartials;

        ConcurrentMap<Region, LikelihoodCore> tmpLikelihoodCores = likelihoodCores;
        likelihoodCores = storedLikelihoodCores;
        storedLikelihoodCores = tmpLikelihoodCores;

        ConcurrentMap<Region, List<Integer>> tmpConstantPatterns = constantPatterns;
        constantPatterns = storedConstantPatterns;
        storedConstantPatterns = tmpConstantPatterns;

        ConcurrentMap<Region, Double> tmpRegionLogLikelihoods = regionLogLikelihoods;
        regionLogLikelihoods = storedRegionLogLikelihoods;
        storedRegionLogLikelihoods = tmpRegionLogLikelihoods;

        super.restore();
    }
}

