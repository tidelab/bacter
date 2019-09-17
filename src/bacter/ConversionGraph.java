/*
 * Copyright (C) 2013 Tim Vaughan <tgvaughan@gmail.com>
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
package bacter;

import bacter.util.parsers.ExtendedNewickBaseVisitor;
import bacter.util.parsers.ExtendedNewickLexer;
import bacter.util.parsers.ExtendedNewickParser;
import beast.core.*;
import beast.evolution.tree.Node;
import beast.evolution.tree.Tree;
import beast.util.TreeParser;
import org.antlr.v4.runtime.CharStream;
import org.antlr.v4.runtime.CharStreams;
import org.antlr.v4.runtime.CommonTokenStream;
import org.antlr.v4.runtime.tree.ParseTree;

import java.io.PrintStream;
import java.util.*;
import java.util.logging.Level;
import java.util.logging.Logger;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import java.util.stream.Collectors;

/**
 * @author Tim Vaughan <tgvaughan@gmail.com>
 */
@Description("Conversion graph based around the clonal frame.")
@Citation("TG Vaughan, D Welch, AJ Drummond, PJ Biggs, T George, NP French\n" +
        "Inferring Ancestral Recombination Graphs from Bacterial Genomic Data,\n" +
        "Genetics, 205(2):857 (2017)")
public class ConversionGraph extends Tree {
    
    /**
     * Unlike Trees, Conversion graphs require an alignment (or at least
 the length of an alignment) to be specified so that the regions of
 the alignment affected by recombination events can be recorded.
     */
    public Input<List<Locus>> lociInput = new Input<>(
            "locus",
            "Locus associated with graph.",
            new ArrayList<>());

    public Input<String> fromStringInput = new Input<>(
            "fromString",
            "Initialise ARG from string representation.");

    public Input<String> fromExtNewickInput = new Input<>(
            "extendedNewick",
            "Initialise ARG from extended Newick representation.");

    public Input<Boolean> wholeLocusConversionsInput = new Input<>(
            "wholeLocusConversionsOnly",
            "Force region boundaries to coincide with locus boundaries.", false);

    /**
     * List of recombinations on graph.
     */
    protected Map<Locus, List<Conversion>> convs;
    protected Map<Locus, List<Conversion>> storedConvs;
    protected final static List<Conversion> emptyConvList = new ArrayList<>();

    /**
     * Event and region lists.
     */
    protected Map<Locus, RegionList> regionLists;
    protected CFEventList cfEventList;

    protected List<Locus> loci, convertibleLoci;
    protected int totalConvertibleSequenceLength;

    @Override
    public void initAndValidate() {

        convs = new HashMap<>();
        storedConvs = new HashMap<>();

        if (lociInput.get().isEmpty())
                throw new RuntimeException("Must specify at least one locus " +
                        "as an input to ConversionGraph.");

        loci = lociInput.get();

        // Sort alignment list lexographically in order of BEASTObject IDs
        loci.sort(Comparator.comparing(BEASTObject::getID));
        convertibleLoci = loci.stream().filter(Locus::conversionsAllowed).collect(Collectors.toList());

        totalConvertibleSequenceLength = 0;
        for (Locus locus : convertibleLoci) {
            convs.put(locus, new ArrayList<>());
            storedConvs.put(locus, new ArrayList<>());
            totalConvertibleSequenceLength += locus.getSiteCount();
        }
        
        if (fromStringInput.get() != null) {
            fromStringOld(fromStringInput.get());
        }

        if (fromExtNewickInput.get() != null) {
            fromExtendedNewick(fromExtNewickInput.get());
        }

        regionLists = new HashMap<>();
        for (Locus locus : loci)
            regionLists.put(locus, new RegionList(this, locus));

        cfEventList = new CFEventList(this);

        super.initAndValidate();
    }

    /**
     * Retrieve total length of all sequence alignments.
     *
     * @return total sequence length
     */
    public int getTotalConvertibleSequenceLength() {
        return totalConvertibleSequenceLength;
    }

    /**
     * Retrieve locus associated with this ACG having a given BEASTObject ID.
     *
     * @param id id of alignment to retrieve
     * @return locus, or null if no locus matches.
     */
    public Locus getLocusByID(String id) {
        for (Locus locus : loci)
            if (locus.getID().equals(id))
                return locus;

        return null;
    }

    /**
     * Retrieve alignments associated with this conversion graph.
     *
     * @return list of associated alignments
     */
    public List<Locus> getConvertibleLoci() {
        return convertibleLoci;
    }

    public boolean wholeLocusModeOn() {
        return wholeLocusConversionsInput.get();
    }

    public void setWholeLocusMode(boolean newVal) {
        wholeLocusConversionsInput.setValue(newVal, this);
    }
    
    /**
     * Add conversion to graph, ensuring conversion list
     * remains sorted.
     *
     * @param conv conversion to add
     */
    public void addConversion(Conversion conv) {
        startEditing(null);
        
        conv.setConversionGraph(this);

        Locus locus = conv.getLocus();

        if (!locus.conversionsAllowed())
            throw new IllegalStateException("Tried to add a conversion to a conversion-free locus.");

        int i;
        for (i=0; i<convs.get(locus).size(); i++)
            if (convs.get(locus).get(i).startSite>conv.startSite)
                break;
        
        convs.get(locus).add(i, conv);
    }
    
    /**
     * Remove recombination from graph.
     *
     * @param conv conversion to remove.
     */
    public void deleteConversion(Conversion conv) {
        startEditing(null);

        if (!conv.getLocus().conversionsAllowed())
            throw new IllegalStateException("Tried to remove a conversion from a conversion-free locus.");

        convs.get(conv.getLocus()).remove(conv);
    }
    
    /**
     * Retrieve list of conversions associated with given locus.
     *
     * @param locus locus with which conversions are associated
     * @return List of conversions.
     */
    public List<Conversion> getConversions(Locus locus) {
        if (locus.conversionsAllowed())
            return convs.get(locus);
        else
            return emptyConvList;
    }

    /**
     * Obtain number of conversion events associated with given locus.
     *
     * @param locus locus with which conversions are associated.
     * @return Number of conversions.
     */
    public int getConvCount(Locus locus) {
        if (locus.conversionsAllowed())
            return convs.get(locus).size();
        else
            return 0;
    }

    /**
     * Obtain total number of conversion events.
     *
     * @return Number of conversions.
     */
    public int getTotalConvCount() {
        int convCount = 0;
        for (Locus locus : convertibleLoci)
            convCount += convs.get(locus).size();

        return convCount;
    }

    /**
     * Obtain index of conversion when conversions are listed in order
     * of alignment and start site.
     *
     * @param conv conversion whose index is required
     * @return Conversion index
     */
    public int getConversionIndex(Conversion conv) {
        int index = 0;
        for (Locus locus : getConvertibleLoci()) {
            if (locus == conv.getLocus()) {
                index += getConversions(locus).indexOf(conv);
                break;
            } else
                index += getConvCount(locus);
        }

        return index;
    }

    /**
     * Get list of contiguous regions having fixed marginal trees
     * associated with given locus.
     *
     * @param locus locus with which regions are associated
     * @return list of regions
     */
    public List<Region> getRegions(Locus locus) {
        return regionLists.get(locus).getRegions();
    }

    /**
     * Obtain number of contiguous single-tree regions associated with
     * given locus.
     *
     * @param locus locus with which regions are associated
     * @return Number of regions.
     */
    public int getRegionCount(Locus locus) {
        return regionLists.get(locus).getRegions().size();
    }

    public int getUselessConvCount() {
        AffectedSiteList asList = new AffectedSiteList(this);

        int count = 0;
        for (int asCount : asList.affectedSiteCount.values()) {
            if (asCount == 0)
                count += 1;
        }

        return count;
    }

    /**
     * Obtain ordered list of events that make up the clonal frame.  Used
     * for ACG probability density calculations and for various state proposal
     * operators.
     * 
     * @return List of events.
     */
    public List<CFEventList.Event> getCFEvents() {
        return cfEventList.getCFEvents();
    }

    /**
     * @return Total length of all edges in clonal frame.
     */
    public double getClonalFrameLength() {
        double length = 0.0;
        for (Node node : m_nodes) {
            if (node.isRoot())
                continue;
            length += node.getLength();
        }
        
        return length;
    }

    /**
     * Check validity of conversions.  Useful for probability densities
     * over the ACG to decide whether to return 0 based on an unphysical
     * state.
     * 
     * @return true if all conversions are valid w.r.t. clonal frame.
     */
    public boolean isInvalid() {
        for (Locus locus : convertibleLoci) {
            for (Conversion conv : convs.get(locus)) {
                if (!conv.isValid()) {
                    return true;
                }
            }
        }
        
        return false;
    }

    /**
     * Produces an extended Newick representation of this ACG.  This
     * method is also used to serialize the state to a state file.
     *
     * @return an extended Newick representation of ACG.
     */
    @Override
    public String toString() {
        return getExtendedNewick();
    }

    /**
     * Produces a string representing the ACG.  The string is
     * composed of a Newick representation of the CF, with additional
     * annotations describing the conversions.  This method is
     * no longer used for state serialization.
     *
     * @return string representation of the ACG
     */
    @Deprecated
    public String toStringOld() {
        StringBuilder sb = new StringBuilder();

        for (Locus locus : convertibleLoci) {
            for (Conversion conv : getConversions(locus)) {
                sb.append(String.format("[&%s,%d,%d,%s,%d,%d,%s] ",
                        locus.getID(),
                        conv.node1.getNr(),
                        conv.startSite,
                        String.valueOf(conv.height1),
                        conv.node2.getNr(),
                        conv.endSite,
                        String.valueOf(conv.height2)));
            }
        }
        sb.append(super.toString());

        return sb.toString();
    }
    
    /**
     * Load ACG from old string representation.
     *
     * @param str string representation of ACG
     */
    @Deprecated
    public void fromStringOld(String str) {
        
        // Extract clonal frame and recombination components of string
        Pattern cfPattern = Pattern.compile("^[^\\(]*(\\(.*)$");
        Matcher cfMatcher = cfPattern.matcher(str);
        
        if (!cfMatcher.find())
            throw new RuntimeException("Error parsing ACG state string.");
        
        // Process clonal frame
        String sNewick = cfMatcher.group(cfMatcher.groupCount());
        try {
            TreeParser parser = new TreeParser();
            parser.thresholdInput.setValue(1e-10, parser);
            parser.offsetInput.setValue(0, parser);
            setRoot(parser.parseNewick(sNewick));
        } catch (Exception ex) {
            Logger.getLogger(ConversionGraph.class.getName()).log(Level.SEVERE, null, ex);
        }

        initArrays();
        
        Pattern convPattern = Pattern.compile("\\[&([^]]*)]");
        Matcher convMatcher = convPattern.matcher(str);
        
        // Process recombinations
        for (Locus locus : convertibleLoci)
            convs.get(locus).clear();

        while(convMatcher.find()) {
            String [] elements = convMatcher.group(1).split(",");

            Locus locus = getLocusByID(elements[0]);
            if (locus == null)
                throw new RuntimeException("Uknown locus id "
                        + elements[0] + ".  Aborting.");

            Node node1 = getNode(Integer.parseInt(elements[1]));
            int startLocus = Integer.parseInt(elements[2]);
            double height1 = Double.parseDouble(elements[3]);
            
            Node node2 = getNode(Integer.parseInt(elements[4]));
            int endLocus = Integer.parseInt(elements[5]);
            double height2 = Double.parseDouble(elements[6]);

            Conversion conv = new Conversion(
                    node1, height1,
                    node2, height2,
                    startLocus, endLocus, this, locus);
            
            addConversion(conv);
        }

        if (isInvalid()) {
            throw new IllegalArgumentException(
                    "Invalid ACG read from string. Aborting.");
        }
    }
    
    @Override
    public void fromXML(final org.w3c.dom.Node node) {
        fromExtendedNewick(node.getTextContent().replaceAll("&amp", "&"));
    }

    @Override
    public ConversionGraph copy() {
        ConversionGraph acg = new ConversionGraph();

        acg.setID(getID());

        acg.index = index;
        acg.root = root.copy();
        acg.nodeCount = nodeCount;
        acg.internalNodeCount = internalNodeCount;
        acg.leafNodeCount = leafNodeCount;

        acg.initArrays();

        acg.m_taxonset.setValue(m_taxonset.get(), acg);
        
        acg.convs = new HashMap<>();
        acg.storedConvs = new HashMap<>();

        acg.loci = loci;
        acg.convertibleLoci = convertibleLoci;
        for (Locus locus : convertibleLoci) {
            acg.convs.put(locus, new ArrayList<>());
            for (Conversion conv : convs.get(locus)) {
                Conversion convCopy = conv.getCopy();
                convCopy.setConversionGraph(acg);
                convCopy.setNode1(acg.m_nodes[conv.getNode1().getNr()]);
                convCopy.setNode2(acg.m_nodes[conv.getNode2().getNr()]);
                acg.convs.get(locus).add(convCopy);
            }

            acg.storedConvs.put(locus, new ArrayList<>());
            for (Conversion conv : storedConvs.get(locus)) {
                Conversion convCopy = conv.getCopy();
                convCopy.setConversionGraph(acg);
                convCopy.setNode1(acg.m_nodes[conv.getNode1().getNr()]);
                convCopy.setNode2(acg.m_nodes[conv.getNode2().getNr()]);
                acg.storedConvs.get(locus).add(convCopy);
            }
        }

        return acg;
    }

    private void generalAssignFrom(StateNode other, boolean fragileAssignment) {
        if (fragileAssignment)
            super.assignFromFragile(other);
        else
            super.assignFrom(other);

        if (other instanceof ConversionGraph) {
            ConversionGraph acg = (ConversionGraph)other;

            loci = acg.loci;
            convertibleLoci = acg.convertibleLoci;

            convs.clear();
            storedConvs.clear();
            for (Locus locus : convertibleLoci) {
                convs.put(locus, new ArrayList<>());
                storedConvs.put(locus, new ArrayList<>());
                for (Conversion conv : acg.getConversions(locus)) {
                    Conversion convCopy = conv.getCopy();
                    convCopy.setConversionGraph(this);
                    convCopy.setNode1(m_nodes[conv.getNode1().getNr()]);
                    convCopy.setNode2(m_nodes[conv.getNode2().getNr()]);
                    convs.get(locus).add(convCopy);
                }
            }

            if (cfEventList == null)
                cfEventList = new CFEventList(this);

            regionLists.clear();
            for (Locus locus : loci) {
                regionLists.put(locus, new RegionList(this, locus));
            }
        }
    }

    /**
     * Use another StateNode to configure this ACG.  If the other StateNode
     * is merely a tree, only the clonal frame is configured.
     * 
     * @param other StateNode used to configure ACG
     */
    @Override
    public void assignFrom(StateNode other) {
        generalAssignFrom(other, false);
    }
    
    @Override
    public void assignFromFragile(StateNode other) {
        generalAssignFrom(other, true);
    }
    
    /**
     * Obtain extended Newick representation of ACG.  Includes Nexus metadata
     * on hybrid leaf nodes describing the alignment sites affected by the
     * conversion event.
     * 
     * @return Extended Newick string.
     */
    public String getExtendedNewick() {
        return getExtendedNewick(true);
    }

    /**
     * Obtain extended Newick representation of ACG.  Optionally Nexus metadata
     * on hybrid leaf nodes describing the alignment sites affected by the
     * conversion event.
     *
     * @param computeAffectedSites if true, compute affected sites
     * @return Extended Newick string.
     */
    public String getExtendedNewick(boolean computeAffectedSites) {
        AffectedSiteList asList = null;
        if (computeAffectedSites)
            asList = new AffectedSiteList(this);

        return extendedNewickTraverse(root, asList, false) + ";";

    }

    /**
     * Obtain extended Newick representation of ACG, including only those
     * conversions which attach to CF edges above non-root nodes.
     * Includes Nexus metadata on hybrid leaf nodes describing the alignment
     * sites affected by the conversion event.
     *
     * @return Extended Newick string.
     */
    public String getTrimmedExtendedNewick() {

        return extendedNewickTraverse(root, null, true) + ";";
    }
    
    private String extendedNewickTraverse(Node node,
                                          AffectedSiteList asList,
                                          boolean intraCFOnly) {
        StringBuilder sb = new StringBuilder();
        
        // Determine sequence of events along this node.
        class Event {
            boolean isArrival;
            double time;
            Conversion conv;
            
            public Event(boolean isArrival, double time, Conversion conv) {
                this.isArrival = isArrival;
                this.time = time;
                this.conv = conv;
            }
        }
        List<Event> events = new ArrayList<>();
        for (Locus locus : getConvertibleLoci()) {
            for (Conversion conv : getConversions(locus)) {
                if (intraCFOnly && conv.node2.isRoot())
                    continue;

                if (conv.node1 == node)
                    events.add(new Event(false, conv.getHeight1(), conv));
                if (conv.node2 == node)
                    events.add(new Event(true, conv.getHeight2(), conv));
            }
        }
        
        // Sort events from oldest to youngest.
        events.sort((Event e1, Event e2) -> {
            if (e1.time > e2.time)
                return -1;
            else
                return 1;
        });

        // Process events.
        
        int cursor = 0;
        
        double lastTime;
        if (node.isRoot())
            lastTime = Double.POSITIVE_INFINITY;
        else
            lastTime = node.getParent().getHeight();

        for (Event event : events) {

            double thisLength;
            if (Double.isInfinite(lastTime))
                thisLength = 0.0;
            else
                thisLength = lastTime - event.time;

            if (event.isArrival) {
                String meta =  String.format(Locale.ENGLISH,
                        "[&conv=%d, region={%d,%d}, locus=\"%s\", relSize=%g",
                        convs.get(event.conv.getLocus()).indexOf(event.conv),
                        event.conv.getStartSite(),
                        event.conv.getEndSite(),
                        event.conv.getLocus().getID(),
                        event.conv.getSiteCount()/(double)event.conv.getLocus().getSiteCount()
                );

                if (asList != null) {
                    meta += String.format(Locale.ENGLISH,
                            ", affectedSites=%d, uselessSiteFraction=%g",
                            asList.affectedSiteCount.get(event.conv),
                            1.0-asList.affectedSiteFraction.get(event.conv));
                }

                if (event.conv.newickMetaDataMiddle != null)
                    meta += ", " + event.conv.newickMetaDataMiddle;

                meta += "]";

                String parentMeta;
                if (event.conv.newickMetaDataTop != null)
                    parentMeta = "[&" + event.conv.newickMetaDataTop + "]";
                else
                    parentMeta = "";

                sb.insert(cursor, "(,#" + getConversionIndex(event.conv)
                        + meta
                        + ":" + (event.conv.height2-event.conv.height1)
                        + ")"
                        + parentMeta
                        + ":" + thisLength);
                cursor += 1;
            } else {
                String meta;
                if (event.conv.newickMetaDataBottom != null)
                    meta = "[&" + event.conv.newickMetaDataBottom + "]";
                else
                    meta = "";

                sb.insert(cursor, "()#" + getConversionIndex(event.conv)
                        + meta
                        + ":" + thisLength);
                cursor += 1;
            }
            
            lastTime = event.time;
        }
        
        // Process this node and its children.

        if (!node.isLeaf()) {
            String subtree1 = extendedNewickTraverse(node.getChild(0), asList, intraCFOnly);
            String subtree2 = extendedNewickTraverse(node.getChild(1), asList, intraCFOnly);
            sb.insert(cursor, "(" + subtree1 + "," + subtree2 + ")");
            cursor += subtree1.length() + subtree2.length() + 3;
        }

        double thisLength;
        if (Double.isInfinite(lastTime))
            thisLength = 0.0;
        else
            thisLength = lastTime - node.getHeight();
        sb.insert(cursor, (node.getNr() + taxaTranslationOffset)
                + node.getNewickMetaData() + ":" + thisLength);
        
        return sb.toString();
    }

    /**
     * Read in an ACG from a string in extended newick format.  Assumes
     * that the network is stored with exactly the same metadata as written
     * by the getExtendedNewick() method.
     *
     * @param string extended newick representation of ACG
     */
    public void fromExtendedNewick(String string) {
        fromExtendedNewick(string, false, taxaTranslationOffset);
    }
    /**
     * Read in an ACG from a string in extended newick format.  Assumes
     * that the network is stored with exactly the same metadata as written
     * by the getExtendedNewick() method.
     *
     * @param string extended newick representation of ACG
     * @param numbered true indicates that the ACG is numbered.
     */
    public void fromExtendedNewick(String string, boolean numbered, int nodeNumberoffset) {

        // Spin up ANTLR
        CharStream input = CharStreams.fromString(string);
        ExtendedNewickLexer lexer = new ExtendedNewickLexer(input);
        CommonTokenStream tokens = new CommonTokenStream(lexer);
        ExtendedNewickParser parser = new ExtendedNewickParser(tokens);
        ParseTree parseTree = parser.tree();

        Map<String, Conversion> convIDMap = new HashMap<>();
        Node root = new ExtendedNewickBaseVisitor<Node>() {

            /**
             * Convert branch lengths to node heights for all nodes in clade.
             *
             * @param node clade parent
             * @return minimum height assigned in clade.
             */
            private double branchLengthsToHeights(Node node) {
                if (node.isRoot())
                    node.setHeight(0.0);
                else
                    node.setHeight(node.getParent().getHeight() - node.getHeight());

                double minHeight = node.getHeight();

                for (Node child : node.getChildren()) {
                    minHeight = Math.min(minHeight, branchLengthsToHeights(child));
                }

                return minHeight;
            }

            /**
             * Remove height offset from all nodes in clade
             * @param node parent of clade
             * @param offset offset to remove
             */
            private void removeOffset(Node node, double offset) {
                node.setHeight(node.getHeight() - offset);

                for (Node child : node.getChildren())
                    removeOffset(child, offset);
            }

            private Node getTrueNode(Node node) {
                if (node.isLeaf()) {
                    assert !convIDMap.containsKey(node.getID());
                    return node;
                }

                if (convIDMap.containsKey(node.getID()))
                    return getTrueNode(node.getChild(0));

                int hybridIdx = -1;
                int nonHybridIdx = -1;
                for (int i=0; i<node.getChildCount(); i++) {
                    if (node.getChild(i).isLeaf() && convIDMap.containsKey(node.getChild(i).getID()))
                        hybridIdx = i;
                    else
                        nonHybridIdx = i;
                }

                if (hybridIdx>0)
                    return getTrueNode(node.getChild(nonHybridIdx));

                return node;
            }

            /**
             * Traverse the newly constructed tree looking for
             * hybrid nodes and using these to set the heights of
             * Conversion objects.
             *
             * @param node parent of clade
             */
            private void findConversionAttachments(Node node) {
                if (convIDMap.containsKey(node.getID())) {
                    Conversion conv = convIDMap.get(node.getID());
                    if (node.isLeaf()) {
                        conv.setHeight1(node.getHeight());
                        conv.setHeight2(node.getParent().getHeight());
                        conv.setNode2(getTrueNode(node.getParent()));
                    } else
                        conv.setNode1(getTrueNode(node));
                }

                for (Node child : node.getChildren())
                    findConversionAttachments(child);
            }

            /**
             * Remove all conversion-associated nodes, leaving only
             * the clonal frame.
             *
             * @param node parent of clade
             * @return new parent of same clade
             */
            private Node stripHybridNodes(Node node) {
                Node trueNode = getTrueNode(node);
                List<Node> trueChildren = new ArrayList<>();

                for (Node child : trueNode.getChildren()) {
                    trueChildren.add(stripHybridNodes(child));
                }

                trueNode.removeAllChildren(false);
                for (Node trueChild : trueChildren)
                    trueNode.addChild(trueChild);

                return trueNode;
            }

            private int numberInternalNodes(Node node, int nextNr) {
                if (node.isLeaf())
                    return nextNr;

                for (Node child : node.getChildren())
                    nextNr = numberInternalNodes(child, nextNr);

                node.setNr(nextNr);

                return nextNr + 1;
            }


            @Override
            public Node visitTree(ExtendedNewickParser.TreeContext ctx) {
                Node root =  visitNode(ctx.node());

                double minHeight = branchLengthsToHeights(root);
                removeOffset(root, minHeight);

                findConversionAttachments(root);

                root = stripHybridNodes(root);
                root.setParent(null);

                if (!numbered)
                    numberInternalNodes(root, root.getAllLeafNodes().size());

                return root;
            }

            @Override
            public Node visitNode(ExtendedNewickParser.NodeContext ctx) {
                Node node = new Node();

                if (ctx.post().hybrid() != null) {
                    String convID = ctx.post().hybrid().getText();
                    node.setID(convID);

                    Conversion conv;
                    if (convIDMap.containsKey(convID))
                        conv = convIDMap.get(convID);
                    else {
                        conv = new Conversion();
                        convIDMap.put(convID, conv);
                    }

                    if (ctx.node().isEmpty()) {
                        String locusID;
                        for (ExtendedNewickParser.AttribContext attribCtx : ctx.post().meta().attrib()) {
                            switch (attribCtx.attribKey.getText()) {
                                case "region":
                                    conv.setStartSite(Integer.parseInt(
                                            attribCtx.attribValue().vector().attribValue(0).getText()));
                                    conv.setEndSite(Integer.parseInt(
                                            attribCtx.attribValue().vector().attribValue(1).getText()));
                                    break;

                                case "locus":
                                    locusID = attribCtx.attribValue().getText();
                                    if (locusID.startsWith("\""))
                                        locusID = locusID.substring(1,locusID.length()-1);

                                    Locus locus = null;
                                    for (Locus thisLocus : getConvertibleLoci()) {
                                        if (thisLocus.getID().equals(locusID))
                                            locus = thisLocus;
                                    }

                                    if (locus == null)
                                        throw new IllegalArgumentException(
                                                "Locus with ID " + locusID + " not found.");

                                    conv.setLocus(locus);
                                    break;

                                default:
                                    break;
                            }
                        }
                    }
                }

                for (ExtendedNewickParser.NodeContext childCtx : ctx.node())
                    node.addChild(visitNode(childCtx));

                if (ctx.post().label() != null) {
                    node.setID(ctx.post().label().getText());
                    node.setNr(Integer.parseInt(ctx.post().label().getText())
                            - nodeNumberoffset);
                }

                node.setHeight(Double.parseDouble(ctx.post().length.getText()));

                return node;
            }
        }.visit(parseTree);

        m_nodes = root.getAllChildNodesAndSelf().toArray(m_nodes);
        nodeCount = m_nodes.length;
        leafNodeCount = root.getAllLeafNodes().size();

        setRoot(root);
        initArrays();

        for (Locus locus : getConvertibleLoci())
            convs.get(locus).clear();

        for (Conversion conv : convIDMap.values())
            addConversion(conv);

    }

    /*
    * StateNode implementation
    */
    
    @Override
    protected void store () {
        super.store();
        
        for (Locus locus : convertibleLoci) {
            storedConvs.get(locus).clear();

            for (Conversion conv : convs.get(locus)) {
                Conversion convCopy = new Conversion();

                convCopy.setLocus(conv.getLocus());
                convCopy.setStartSite(conv.getStartSite());
                convCopy.setEndSite(conv.getEndSite());
                convCopy.setHeight1(conv.getHeight1());
                convCopy.setHeight2(conv.getHeight2());
                convCopy.newickMetaDataBottom = conv.newickMetaDataBottom;
                convCopy.newickMetaDataMiddle = conv.newickMetaDataMiddle;
                convCopy.newickMetaDataTop = conv.newickMetaDataTop;

                convCopy.setNode1(m_storedNodes[conv.getNode1().getNr()]);
                convCopy.setNode2(m_storedNodes[conv.getNode2().getNr()]);

                convCopy.setConversionGraph(this);

                storedConvs.get(locus).add(convCopy);
            }
        }
    }
    
    @Override
    public void restore() {
        super.restore();
        
        Map<Locus, List<Conversion>> tmp = storedConvs;
        storedConvs = convs;
        convs = tmp;

        cfEventList.makeDirty();
        for (Locus locus : loci)
            regionLists.get(locus).makeDirty();
    }

    @Override
    public void startEditing(Operator operator) {
        if (state != null)
            super.startEditing(operator);

        if (cfEventList != null)
            cfEventList.makeDirty();

        if (regionLists != null)
            for (RegionList regionList : regionLists.values())
                regionList.makeDirty();
    }

    /**
     * @return true iff clonal frame is dirty
     */
    public boolean clonalFrameIsDirty() {
        for (Node node : getNodesAsArray())
            if (node.isDirty() > Tree.IS_CLEAN)
                return true;

        return false;
    }
    
    /*
     * Loggable implementation.
     */
    @Override
    public void init(PrintStream out) {
        Node node = getRoot();
        out.println("#NEXUS\n");
        out.println("Begin taxa;");
        out.println("\tDimensions ntax=" + getLeafNodeCount() + ";");
        out.println("\t\tTaxlabels");
        printTaxa(node, out, getNodeCount() / 2);
        out.println("\t\t\t;");
        out.println("End;\n");

        out.println("Begin bacter;");
        out.print("\tloci");
        for (Locus locus : loci)
            out.print(" " + locus.getID() + ":" + locus.getSiteCount());
        out.println(";\nEnd;\n");

        out.println("Begin trees;");
        out.println("\tTranslate");
        printTranslate(node, out, getNodeCount() / 2);
        out.print(";");
    }

    @Override
    public void log(long nSample, PrintStream out) {
        ConversionGraph arg = (ConversionGraph) getCurrent();
        
        out.print(String.format("tree STATE_%d = [&R] %s",
                nSample, arg.getExtendedNewick()));
    }
}
