/*
 * Copyright (C) 2016 Tim Vaughan <tgvaughan@gmail.com>
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

package bacter.devutils;

import bacter.Conversion;
import bacter.ConversionGraph;
import beast.evolution.tree.TreeUtils;
import bacter.Locus;
import bacter.acgannotator.ACGAnnotator;
import bacter.model.ACGLikelihood;
import bacter.util.ACGLogReader;
import bacter.util.BacterACGLogReader;
import bacter.util.COACGLogFileReader;
import beast.evolution.alignment.Alignment;
import beast.evolution.datatype.DataType;
import beast.evolution.sitemodel.SiteModel;
import beast.evolution.substitutionmodel.EigenDecomposition;
import beast.evolution.substitutionmodel.Frequencies;
import beast.evolution.substitutionmodel.HKY;
import beast.evolution.substitutionmodel.SubstitutionModel;
import beast.evolution.tree.Node;
import beast.util.NexusParser;

import javax.xml.stream.XMLStreamException;
import java.io.File;
import java.io.IOException;
import java.io.PrintStream;
import java.util.*;

/**
 * @author Tim Vaughan <tgvaughan@gmail.com>
 */
public class DifferenceFromTrueACG {

    private static class Options {
        double burninPerc = 10.0;
        double boundaryTol = 0.25;
        double ageTol = 0.25;
        double minSupportClade = 0.5;
        double minSupportConv = 0.5;
        File logFile, truthFile, aliFile, summarizedAcgOut, outFileSummarizedClade, outFileSummarizedConv, outFileTrueConv;
        boolean useCOFormat = false;
    }

    public static void printUsageAndExit(int exitCode) {
        System.out.println("Usage: DifferenceFromTrueACG [-burnin b] [-boundaryTol t] [-ageTol t] [-minSupportClade t] [-minSupportConv t] [-co] truth.tree log.trees true_alignment_file output_summarized_ACG output_compare_summarized_clades output_compare_summarized_conversions output_compare_true_conversions");
        System.exit(exitCode);
    }

    /**
     * Process command line arguments.
     *
     * @param args
     * @return
     */
    public static Options processArguments(String[] args) {

        Options options = new Options();

        int i = 0;
        while (i < args.length && args[i].startsWith("-")) {
            switch (args[i].substring(1)) {
                case "co":
                    options.useCOFormat = true;
                    break;

                case "burnin":
                    i += 1;
                    if (i >= args.length)
                        printUsageAndExit(1);
                    try {
                        options.burninPerc = Double.valueOf(args[i]);
                    } catch (NumberFormatException e) {
                        System.out.println("Argument to -burnin must be a number.");
                        printUsageAndExit(1);
                    }
                    break;

                case "boundaryTol":
                    i += 1;
                    if (i >= args.length)
                        printUsageAndExit(1);
                    try {
                        options.boundaryTol = Double.valueOf(args[i]);
                    } catch (NumberFormatException e) {
                        System.out.println("Argument to -boundaryTol must be a number.");
                        printUsageAndExit(1);
                    }
                    break;

                case "minSupportClade":
                    i += 1;
                    if (i >= args.length)
                        printUsageAndExit(1);
                    try {
                        options.minSupportClade = Double.valueOf(args[i]);
                    } catch (NumberFormatException e) {
                        System.out.println("Argument to -minSupportClade must be a number.");
                        printUsageAndExit(1);
                    }
                    break;

                case "minSupportConv":
                    i += 1;
                    if (i >= args.length)
                        printUsageAndExit(1);
                    try {
                        options.minSupportConv = Double.valueOf(args[i]);
                    } catch (NumberFormatException e) {
                        System.out.println("Argument to -minSupportConv must be a number.");
                        printUsageAndExit(1);
                    }
                    break;

                case "ageTol":
                    i += 1;
                    if (i >= args.length)
                        printUsageAndExit(1);
                    try {
                        options.ageTol = Double.valueOf(args[i]);
                    } catch (NumberFormatException e) {
                        System.out.println("Argument to -ageTol must be a number.");
                        printUsageAndExit(1);
                    }
                    break;

                default:
                    System.err.println("Unknown argument: " + args[i]);
                    printUsageAndExit(1);
            }

            i++;
        }

        if (args.length - i < 7)
            printUsageAndExit(0);

        options.truthFile = new File(args[i++]);
        options.logFile = new File(args[i++]);
        options.aliFile = new File(args[i++]);
        options.summarizedAcgOut = new File(args[i++]);
        options.outFileSummarizedClade = new File(args[i++]);
        options.outFileSummarizedConv = new File(args[i++]);
        options.outFileTrueConv = new File(args[i]);

        return options;
    }

    public static class Clade extends BitSet {
        public double age;
    }

    public static Clade getClades(Clade[] clades, Node node) {
        Clade clade = new Clade();
        clade.age = node.getHeight();

        if (node.isLeaf()) {
            clade.set(node.getNr());
        } else {
            for (Node child : node.getChildren())
                clade.or(getClades(clades, child));
        }

        clades[node.getNr()] = clade;

        return clade;
    }

//    /**
//     * Count number of true clades which exist in the provided sampled ARG.
//     *
//     * @param trueClades clades in true arg
//     * @param clades clades in sampled arg
//     * @param ageTol maximum relative age error
//     * @param cladeHist
//     * @return number of found clades
//     */
//    public static int countSampledTrueClades(Clade[] trueClades, Clade[] clades, double ageTol,
//                                             Map<Clade, Integer> cladeHist) {
//        int foundClades = 0;
//        for (Clade trueClade : trueClades) {
//            for (Clade clade : clades) {
//                if (!clade.equals(trueClade)) // || (clade.cardinality()>1 && (trueClade.age - clade.age)/trueClade.age > ageTol))
//                    continue;
//
//                foundClades += 1;
//
//                cladeHist.put(clade, cladeHist.get(clade)+1);
//
//                break;
//            }
//        }
//
//        return foundClades;
//    }
//
//    /**
//     * Count the number of true conversions which have correspondences on the
//     * provided sampled ARG.
//     *
//     * @param trueACG true arg
//     * @param trueClades clades in true arg
//     * @param acg sampled arg
//     * @param clades clades in sampled arg
//     * @param boundaryTol minimum relative error in region boundaries to allow.
//     * @param convHist
//     * @return number of found conversions
//     */
//    public static int countSampledTrueConversions(ConversionGraph trueACG, Clade[] trueClades,
//                                                  ConversionGraph acg, Clade[] clades,
//                                                  double boundaryTol, double ageTol,
//                                                  Map<Conversion, Integer> convHist) {
//        int count = 0;
//        for (Conversion trueConv : trueACG.getConversions(trueACG.getConvertibleLoci().get(0))) {
//            Clade trueFromClade = trueClades[trueConv.getNode1().getNr()];
//            Clade trueToClade = trueClades[trueConv.getNode2().getNr()];
//            for (Conversion conv : acg.getConversions(acg.getConvertibleLoci().get(0)))  {
//                Clade fromClade = clades[conv.getNode1().getNr()];
//                Clade toClade = clades[conv.getNode2().getNr()];
//
//                if (fromClade.equals(trueFromClade) && toClade.equals(trueToClade)) {
//                    if (    1.0*Math.abs(trueConv.getStartSite()-conv.getStartSite())/trueConv.getSiteCount() <= boundaryTol
//                            && 1.0*Math.abs(trueConv.getEndSite()-conv.getEndSite())/trueConv.getSiteCount() <= boundaryTol
////                            && Math.abs(trueConv.getHeight1()-conv.getHeight1())/trueConv.getHeight1() <= ageTol
//                            //&& Math.abs(trueConv.getHeight2()-conv.getHeight2())/trueConv.getHeight2() <= ageTol
//                            )
//                    {
//                        count += 1;
//
//                        convHist.put(trueConv, convHist.get(trueConv) + 1);
//
//                        break;
//                    }
//                }
//            }
//        }
//
//        return count;
//    }
//
//    /**
//     * Count the number of actual clades present with given support in posterior.
//     *
//     * @param cladeHist Mapping from clades to the number of sampled ACGs that
//     *                  contain them.
//     * @param nSampledACGs Total number of sampled ACGs
//     * @param minimumSupport Minimum fraction of samples clades must appear to
//     *                       be considered "recovered"
//     * @return number of recovered clades
//     */
//    public static int countRecoveredTrueClades(Map<Clade, Integer> cladeHist,
//                                               int nSampledACGs, double minimumSupport) {
//        int count = 0;
//
//        int sampleThresh = (int)Math.round(nSampledACGs*minimumSupport);
//
//        for (Clade clade : cladeHist.keySet()) {
//            if (cladeHist.get(clade) >= sampleThresh)
//                count += 1;
//        }
//
//        return count;
//    }
//
//    public static int countRecoveredTrueConvs(Map<Conversion, Integer> convHist,
//                                              int nSampledACGs, double minimumSupport) {
//        int count = 0;
//
//        int sampleThresh = (int)Math.round(nSampledACGs*minimumSupport);
//
//        for (Conversion conv : convHist.keySet()) {
//            if (convHist.get(conv) >= sampleThresh)
//                count += 1;
//        }
//
//        return count;
//    }
//
//    // Count the number of clades and conversions recovered in the summarized ACG
//    public static int countRecoveredTrueCladesInSumACG(Clade[] trueClades,
//                                               ConversionGraph summarizedACG) {
//        int foundClades = 0;
//        Clade[] clades = new Clade[summarizedACG.getNodeCount()];
//        getClades(clades, summarizedACG.getRoot());
//        for (Clade trueClade : trueClades) {
//            for (Clade clade : clades) {
//                if (!clade.equals(trueClade)) // || (clade.cardinality()>1 && (trueClade.age - clade.age)/trueClade.age > ageTol))
//                    continue;
//
//                foundClades += 1;
//
//                break;
//            }
//        }
//
//        return foundClades;
//    }

    // Get and print info on summarized ACG clades
    public static void compareSumACGClades(Clade[] trueClades,
                                           ConversionGraph summarizedACG, PrintStream ps) {
        Clade[] clades = new Clade[summarizedACG.getNodeCount()];
        getClades(clades, summarizedACG.getRoot());
        ps.println("CladeNb" + "\t" + "isTrue" + "\t" + "support");
        for (int cladeNum = 0; cladeNum < clades.length; cladeNum++) {
            int isTrue = 0;
            double support = Double.parseDouble(summarizedACG.getNode(cladeNum).metaDataString.replaceAll(".*posterior=([^,]*),.*", "$1"));
            for (Clade trueClade : trueClades) {
                if (!clades[cladeNum].equals(trueClade)) // || (clade.cardinality()>1 && (trueClade.age - clade.age)/trueClade.age > ageTol))
                    continue;
                isTrue = 1;
                break;
            }
            ps.println(cladeNum + "\t" + isTrue + "\t" + support);
        }

    }

    // Get and print info on summarized ACG conversions
    public static void compareSumACGConv(ConversionGraph trueACG, Clade[] trueClades,
                                         ConversionGraph summarizedACG, PrintStream ps) {

        Clade[] clades = new Clade[summarizedACG.getNodeCount()];
        getClades(clades, summarizedACG.getRoot());
        int convNum = 1;
        ps.println("ConvNb" + "\t" + "isTrue" + "\t" + "support" + "\t" + "edgeError");

        for (Conversion conv : summarizedACG.getConversions(summarizedACG.getConvertibleLoci().get(0))) {

            int isTrue = 0;
            double edgeError = .0;

            //get start and end clades of the conversion
            Clade fromClade = clades[conv.getNode1().getNr()];
            Clade toClade = clades[conv.getNode2().getNr()];

            //get start site estimate and 95% HPD bounds
            HashMap<String, Integer> startSite = new HashMap<>();
            HashMap<String, Integer> endSite = new HashMap<>();

            startSite.put("estimate", conv.getStartSite());
            startSite.put("lowHPD", Integer.parseInt(conv.newickMetaDataMiddle.replaceAll(".*startSite_95%_HPD=\\{([0-9]*).*", "$1")));
            startSite.put("upHPD", Integer.parseInt(conv.newickMetaDataMiddle.replaceAll(".*startSite_95%_HPD=\\{[0-9]*,([0-9]*).*", "$1")));

            endSite.put("estimate", conv.getEndSite());
            endSite.put("lowHPD", Integer.parseInt(conv.newickMetaDataMiddle.replaceAll(".*endSite_95%_HPD=\\{([0-9]*).*", "$1")));
            endSite.put("upHPD", Integer.parseInt(conv.newickMetaDataMiddle.replaceAll(".*endSite_95%_HPD=\\{[0-9]*,([0-9]*).*", "$1")));

            //get posterior support
            Double posterior = Double.parseDouble(conv.newickMetaDataMiddle.replaceAll(".*posterior=([^,]*).*", "$1"));

            for (Conversion trueConv : trueACG.getConversions(trueACG.getConvertibleLoci().get(0))) {
                Clade trueFromClade = trueClades[trueConv.getNode1().getNr()];
                Clade trueToClade = trueClades[trueConv.getNode2().getNr()];
                int trueStartSite = trueConv.getStartSite();
                int trueEndSite = trueConv.getEndSite();

                if (fromClade.equals(trueFromClade) && toClade.equals(trueToClade)) {
                    //if (trueStartSite >= startSite.get("lowHPD") && trueStartSite <= startSite.get("upHPD") && trueEndSite >= endSite.get("lowHPD") && trueEndSite <= endSite.get("upHPD")) {
                    //if (1.0*Math.abs(trueConv.getStartSite()-conv.getStartSite())/trueConv.getSiteCount() <= 0.5 && 1.0*Math.abs(trueConv.getEndSite()-conv.getEndSite())/trueConv.getSiteCount() <= 0.5){
                    if (conv.getStartSite() < trueEndSite && conv.getEndSite() > trueStartSite) {
                        isTrue = 1;
                        edgeError = (double) ((Math.abs(trueStartSite - startSite.get("estimate")) + Math.abs(trueEndSite - endSite.get("estimate")))) / (2 * trueConv.getSiteCount());
                        break;
                    }
                }
            }
            ps.println(convNum++ + "\t" + isTrue + "\t" + posterior + "\t" + edgeError);
        }
    }

    //Get and print length of true conversions with their support in the summarized ACG (this is only useful if some conversions were not sampled at all, other wise we could obtain the same info with the previous method

    public static void compareTrueConv(ConversionGraph trueACG, Clade[] trueClades,
                                       ConversionGraph summarizedACG, Alignment trueAlignment, PrintStream ps) {

        //get likelihood of simulated data given true parameter value
        Frequencies freqs = new Frequencies();
        freqs.initByName("data",trueAlignment);
        HKY hkySubstModel = new HKY();
        hkySubstModel.initByName("kappa","4","frequencies",freqs);
        SiteModel siteModel = new SiteModel();
        siteModel.initByName("gammaCategoryCount", "4", "shape","0.2","mutationRate","5e-6","proportionInvariant","0","substModel",hkySubstModel);
        ACGLikelihood acgLikelihoodWithConv = new ACGLikelihood();
        acgLikelihoodWithConv.initByName(
                "locus", trueACG.getLocusByID("locus"),
                "data", trueAlignment,
                "tree", trueACG,
                "siteModel", siteModel);
        double LikelihoodWithConv = acgLikelihoodWithConv.calculateLogP();
        //get clades in summarized ACG
        Clade[] clades = new Clade[summarizedACG.getNodeCount()];
        getClades(clades, summarizedACG.getRoot());

        int convNum = 1;
        ps.println("TrueConvNb" + "\t" + "supportInSummarizedACG" + "\t" + "edgeError" + "\t" + "tractLength" + "\t" + "startHeight" + "\t" + "endHeigth" + "\t" + "likelihoodGain" + "\t" + "isUseless" + "\t" + "distanceBetweenConv");

        for (Conversion trueConv : trueACG.getConversions(trueACG.getConvertibleLoci().get(0))) {
            boolean isUseless = false;
            Clade trueFromClade = trueClades[trueConv.getNode1().getNr()];
            Clade trueToClade = trueClades[trueConv.getNode2().getNr()];
            int trueStartSite = trueConv.getStartSite();
            int trueEndSite = trueConv.getEndSite();
            double edgeError = 1;
            double posterior = 0.0;
            int tractLength = trueConv.getSiteCount();
            double startHeight = trueConv.getHeight2();
            double endHeight = trueConv.getHeight1();
            //get common ancestor of parental nodes
            double convDist;
            Set<String> leafnodes = new HashSet<>();
            if (trueConv.getNode1().isLeaf()){
                leafnodes.add(trueConv.getNode1().getID());
            } else{
                for(Node node:trueConv.getNode1().getAllLeafNodes())
                    leafnodes.add(node.getID());
            }
            if (trueConv.getNode2().isLeaf()){
                leafnodes.add(trueConv.getNode2().getID());
            } else{
                for(Node node:trueConv.getNode2().getAllLeafNodes())
                    leafnodes.add(node.getID());
            }
            Node commonAncestor = TreeUtils.getCommonAncestorNode(trueACG,leafnodes);
            //get evolutionary distance between the departure and the arrival of the conversion
            if (trueConv.getNode1().equals(trueConv.getNode2())){
                isUseless = true;
                convDist = trueConv.getHeight2() - trueConv.getHeight1();
            } else{
                convDist = Math.abs(commonAncestor.getHeight()-trueConv.getHeight1()) + Math.abs(commonAncestor.getHeight()-trueConv.getHeight2());
            }
            //get the gain in tree likelihood brought by this conversion (ratio of tree likelihood with our without the conversion)
            ConversionGraph trueACGWithoutConv = trueACG.copy();
            trueACGWithoutConv.getConversions(trueACGWithoutConv.getLocusByID("locus")).remove(convNum-1);
            ACGLikelihood acgLikelihoodWithoutConv = new ACGLikelihood();
            acgLikelihoodWithoutConv.initByName(
                    "locus", trueACGWithoutConv.getLocusByID("locus"),
                    "data", trueAlignment,
                    "tree", trueACGWithoutConv,
                    "siteModel", siteModel);

            double likelihoodGain=LikelihoodWithConv - acgLikelihoodWithoutConv.calculateLogP();

            //for each true conversion, check if it is present in the summarized ACG and get posterior support as well as tract boundaries relative error
            for (Conversion conv : summarizedACG.getConversions(summarizedACG.getConvertibleLoci().get(0))) {
                //get start and end clades of the conversion
                Clade fromClade = clades[conv.getNode1().getNr()];
                Clade toClade = clades[conv.getNode2().getNr()];

                //get start site estimate and 95% HPD bounds
                HashMap<String, Integer> startSite = new HashMap<>();
                HashMap<String, Integer> endSite = new HashMap<>();

                startSite.put("estimate", conv.getStartSite());
                startSite.put("lowHPD", Integer.parseInt(conv.newickMetaDataMiddle.replaceAll(".*startSite_95%_HPD=\\{([0-9]*).*", "$1")));
                startSite.put("upHPD", Integer.parseInt(conv.newickMetaDataMiddle.replaceAll(".*startSite_95%_HPD=\\{[0-9]*,([0-9]*).*", "$1")));

                endSite.put("estimate", conv.getEndSite());
                endSite.put("lowHPD", Integer.parseInt(conv.newickMetaDataMiddle.replaceAll(".*endSite_95%_HPD=\\{([0-9]*).*", "$1")));
                endSite.put("upHPD", Integer.parseInt(conv.newickMetaDataMiddle.replaceAll(".*endSite_95%_HPD=\\{[0-9]*,([0-9]*).*", "$1")));

                if (fromClade.equals(trueFromClade) && toClade.equals(trueToClade)) {
                    //if (trueStartSite >= startSite.get("lowHPD") && trueStartSite <= startSite.get("upHPD") && trueEndSite >= endSite.get("lowHPD") && trueEndSite <= endSite.get("upHPD")) {
                    //if (1.0*Math.abs(trueConv.getStartSite()-conv.getStartSite())/trueConv.getSiteCount() <= 0.5 && 1.0*Math.abs(trueConv.getEndSite()-conv.getEndSite())/trueConv.getSiteCount() <= 0.5){
                    if (conv.getStartSite() < trueEndSite && conv.getEndSite() > trueStartSite) {
                        posterior = Double.parseDouble(conv.newickMetaDataMiddle.replaceAll(".*posterior=([^,]*).*", "$1"));
                        edgeError = (double) ((Math.abs(trueStartSite - startSite.get("estimate")) + Math.abs(trueEndSite - endSite.get("estimate")))) / (2 * trueConv.getSiteCount());
                        break;
                    }
                }
            }
            ps.println(convNum++ + "\t" + posterior + "\t" + edgeError + "\t" + tractLength + "\t" + startHeight + "\t" + endHeight + "\t" + likelihoodGain + "\t" + isUseless + "\t" + convDist);
        }
    }

    //main

    public static void main(String[] args) throws IOException, XMLStreamException {

        Options options = processArguments(args);

        // Load true ARG

        BacterACGLogReader truthReader = new BacterACGLogReader(options.truthFile, 0);
        if (truthReader.getACGCount() != 1) {
            System.out.println("Expected exactly 1 ACG in truth file. Found " +
                    truthReader.getACGCount());
            System.exit(1);
        }

        ConversionGraph trueACG = null;
        for (ConversionGraph acg : truthReader)
            trueACG = acg;

        // Determine clades present in truth
        Clade[] trueClades = new Clade[trueACG.getNodeCount()];
        getClades(trueClades, trueACG.getRoot());
        Set<Clade> trueCladeSet = new HashSet<>(Arrays.asList(trueClades));

//         // Set up histograms
//
//        Map<Clade, Integer> cladeHist = new HashMap<>();
//        for (Clade clade : trueClades)
//            cladeHist.put(clade, 0);
//
//        Map<Conversion, Integer> convHist = new HashMap<>();
//        for (Conversion conv : trueACG.getConversions(trueACG.getConvertibleLoci().get(0)))
//            convHist.put(conv, 0);

        // Set up ARG log file reader

        ACGLogReader logReader;
        if (options.useCOFormat) {
            logReader = new COACGLogFileReader(options.logFile, options.burninPerc);
        } else {
            logReader = new BacterACGLogReader(options.logFile, options.burninPerc);
        }

        // Load alignments
        Alignment trueAlignment;
        NexusParser nexusParser = new NexusParser();
        nexusParser.parseFile(options.aliFile);
        trueAlignment = nexusParser.m_alignment;

        // Compute summary ACG
        ACGAnnotator acgAnnotator = new ACGAnnotator(options.logFile, options.summarizedAcgOut, options.burninPerc, 1);
        ConversionGraph summarizedACG = acgAnnotator.getSummarizedACG();

        // Compute and write summary statistics to output file

//        // Count number of recovered clades and conversions in sampled ACGs
//        try (PrintStream ps = new PrintStream(options.outFile)) {
//            ps.println("trueCladeCount sampledTrueCladeCount trueConvCount sampledConvCount sampledTrueConvCount");
//
//            for (ConversionGraph acg : logReader) {
//
//                Clade[] clades = new Clade[acg.getNodeCount()];
//                getClades(clades, acg.getRoot());
//
//                List<Double> timeErrors = new ArrayList<>();
//                int sampledTrueClades = countSampledTrueClades(trueClades, clades, options.ageTol, cladeHist);
//                int sampledTrueConvs = countSampledTrueConversions(trueACG, trueClades, acg, clades,
//                        options.boundaryTol, options.ageTol, convHist);
//
//                ps.println(trueACG.getNodeCount() + "\t" +
//                        sampledTrueClades + "\t" +
//                        trueACG.getConvCount(trueACG.getConvertibleLoci().get(0)) + "\t" +
//                        acg.getConvCount(acg.getConvertibleLoci().get(0)) + "\t" +
//                        sampledTrueConvs);
//            }
//        }
//        // Count number of recovered clades and conversions in sampled ACGs with a support threshold
//        try (PrintStream ps = new PrintStream(options.summaryFile)) {
//            ps.println("trueCladeCount recoveredCladeCount trueConvCount recoveredConvCount");
//
//            int recoveredTrueClades = countRecoveredTrueClades(cladeHist,
//                    logReader.getCorrectedACGCount(), options.minSupportClade);
//
//            int recoveredTrueConvs = countRecoveredTrueConvs(convHist,
//                    logReader.getCorrectedACGCount(), options.minSupportConv);
//
//            ps.println(trueACG.getNodeCount() + "\t" +
//                    recoveredTrueClades + "\t" +
//                    trueACG.getConvCount(trueACG.getConvertibleLoci().get(0)) + "\t" +
//                    recoveredTrueConvs);
//        }

        // Assess clades recovered in summarized ACG
        try (PrintStream ps = new PrintStream(options.outFileSummarizedClade)) {
            compareSumACGClades(trueClades, summarizedACG, ps);
        }
        // Assess conversions recovered in summarized ACG
        try (PrintStream ps = new PrintStream(options.outFileSummarizedConv)) {
            compareSumACGConv(trueACG, trueClades, summarizedACG, ps);

        }
        // Check if true conversions were recovered
        try (PrintStream ps = new PrintStream(options.outFileTrueConv)) {
            compareTrueConv(trueACG, trueClades, summarizedACG, trueAlignment, ps);
        }
    }
}
