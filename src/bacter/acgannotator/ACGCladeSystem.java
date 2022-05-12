/*
 * Copyright (C) 2015 Tim Vaughan <tgvaughan@gmail.com>
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

package bacter.acgannotator;

import bacter.Conversion;
import bacter.ConversionGraph;
import bacter.Locus;
import beast.app.treeannotator.CladeSystem;
import beast.evolution.tree.Node;
import beast.util.Randomizer;

import java.util.*;
import java.util.function.BiFunction;

/**
 * Adds conversion summary tools to CladeSystem.
 *
 * @author Tim Vaughan <tgvaughan@gmail.com>
 */
public class ACGCladeSystem extends CladeSystem {

    protected Map<BitSetPair, Map<Locus, List<Conversion>>> conversionLists = new HashMap<>();
    protected Map<BitSetPair, List<Conversion>> conversionListsTemp = new HashMap<>();
    protected List<Map<BitSet, Map<BitSet, Long>>> geneFlow = new ArrayList<>();
    protected BitSet[] bitSets;

    protected int acgIndex = 1;

    public ACGCladeSystem() { }

    public ACGCladeSystem(ConversionGraph acg) {
        add(acg, true);
    }

    /**
     * Assemble list of bitSets for this ACG.
     */
    public BitSet[] getBitSets(ConversionGraph acg) {

        if (bitSets == null)
            bitSets = new BitSet[acg.getNodeCount()];

        applyToClades(acg.getRoot(), (cladeNode, bits) -> {
            bitSets[cladeNode.getNr()] = bits;
            return null;
        });

        return bitSets;
    }

    /**
     * Add conversions described on provided acg to the internal list
     * for later summary.
     *
     * @param acg conversion graph from which to extract conversions
     * @param receiverBranchMode
     */
    public void collectConversions(ConversionGraph acg, boolean receiverBranchMode) {
        getBitSets(acg);

        Map<BitSet,Map<BitSet,Long>> geneFlowTemp = new HashMap<>();

        // Assemble list of conversions for each pair of clades on each locus
        for (Locus locus : acg.getConvertibleLoci()) {

            conversionListsTemp.clear();
            for (Conversion conv : acg.getConversions(locus))  {
                conv.acgIndex = acgIndex;
                BitSetPair bsPair = new BitSetPair(conv);
                if(receiverBranchMode)
                    bsPair.to = new BitSet(); //set the to bitSet to empty if the receiverBranch mode is used //TODO: receiverBranchMode

                if (!conversionListsTemp.containsKey(bsPair))
                    conversionListsTemp.put(bsPair, new ArrayList<>());

                conversionListsTemp.get(bsPair).add(conv);

                // Record gene flow
                if (!geneFlowTemp.containsKey(bsPair.from))
                    geneFlowTemp.put(bsPair.from, new HashMap<>());

                long oldFlow = 0;
                if (geneFlowTemp.get(bsPair.from).containsKey(bsPair.to))
                    oldFlow = geneFlowTemp.get(bsPair.from).get(bsPair.to);

                geneFlowTemp.get(bsPair.from).put(bsPair.to, oldFlow + conv.getSiteCount());
            }

            // Merge overlapping conversions:
            for (BitSetPair bsPair : conversionListsTemp.keySet()) {
                List<Conversion> merged = mergeOverlappingConvs(
                            conversionListsTemp.get(bsPair));

                if (!conversionLists.containsKey(bsPair))
                    conversionLists.put(bsPair, new HashMap<>());
                if (!conversionLists.get(bsPair).containsKey(locus))
                    conversionLists.get(bsPair).put(locus, new ArrayList<>());

                conversionLists.get(bsPair).get(locus).addAll(merged);

            }
        }

        geneFlow.add(geneFlowTemp);

        acgIndex += 1;
    }

    private List<Conversion> mergeOverlappingConvs(List<Conversion> conversions) {
        List<Conversion> mergedList = new ArrayList<>();//TOREMOVE: create empty list of conversions

        List<Conversion> convOrderedByStart = new ArrayList<>(conversions);//TOREMOVE: create list of conversion from input
        convOrderedByStart.sort((o1, o2) -> o1.getStartSite() - o2.getStartSite());//TOREMOVE: order them by start site

        List<Conversion> convOrderedByEnd = new ArrayList<>(conversions);//TOREMOVE: create list of conversion from input
        convOrderedByEnd.sort((o1, o2) -> o1.getEndSite() - o2.getEndSite());//TOREMOVE: order them by end site


        int nActive = 0;//TOREMOVE: set nb. of active conversion to 0
        Conversion currentMergedConv = null;//TOREMOVE: create a null currentMergedConv that will be added to the list
        int mergedConvCount = 0;//TOREMOVE: conv count of currentMergedConv
        List<Double> mergedConvHeight1 = new ArrayList<>();//we use Lists for the mergedConvHeigths so that we don't have to chose their length (this is actually only useful for the heights2 with the receiverBranchMode //TODO: receiverBranchMode
        List<Double> mergedConvHeight2 = new ArrayList<>();
        List<Node> mergedNodes2 = new ArrayList<>();

        //TODO: check adjustment (circular genome)

        int numFirst = 0;//TOREMOVE: circulargenome adjustment: initialize numFirst (nb. of conversions in the first merged conversion in case it is overlapping)
        boolean firstStep = true;//TOREMOVE: circulargenome adjustment: initialize firstep to true
        int minOverlapStart = Integer.MAX_VALUE;//TOREMOVE: circulargenome adjustment
        int indConv = 0;//TOREMOVE: circulargenome adjustment
        for (int i = 0; i < convOrderedByEnd.size(); i++) {//TOREMOVE: circulargenome adjustment: for i in total nb. of conv
            Conversion conv = convOrderedByStart.get(indConv);//TOREMOVE: circulargenome adjustment: get the next conversion ordered by start (starting from the first
            if (conv.getEndSite() < conv.getStartSite()) {//TOREMOVE: circulargenome adjustment: if it is overlapping we add it to the currentMerged conversion
                nActive += 1;//TOREMOVE: circulargenome adjustment: increment nActive
                mergedConvCount += 1;//TOREMOVE: circulargenome adjustment: increment mergedConvcount
                //mergedConvHeight1 += conv.getHeight1();//TOREMOVE: circulargenome adjustment: increment merged heights
                //mergedConvHeight2 += conv.getHeight2();//TOREMOVE: circulargenome adjustment:
                minOverlapStart = Math.min(minOverlapStart, conv.getStartSite());//TOREMOVE: circulargenome adjustment: set min overlap start to the first overlapping conv encountered
                currentMergedConv = conv.getStartSite() <= minOverlapStart ? conv.getCopy() : currentMergedConv;//TOREMOVE: circulargenome adjustment: if this is the first overlapping conversion encountered, we initialise the currentmerged conv with it
                currentMergedConv.acgIndex = conv.acgIndex;//TOREMOVE: circulargenome adjustment: set acgIndex
                convOrderedByStart.remove(indConv);//TOREMOVE: circulargenome adjustment: remove conv from convorderedBystart
                indConv -= 1;//TOREMOVE: circulargenome adjustment: decrement indConv
            }
            indConv += 1;//TOREMOVE: circulargenome adjustment: increment indConv
        }

        while (!convOrderedByStart.isEmpty() || !convOrderedByEnd.isEmpty()) {//TOREMOVE: while we still have some conversions in the remaining part of the genome

            int nextStart = convOrderedByStart.isEmpty()//TOREMOVE: set nexstart to the start of the first element of convOrderedByStart (or MAX_VALUE IF IT IS EMPTY)
                    ? Integer.MAX_VALUE
                    : convOrderedByStart.get(0).getStartSite();

            int nextEnd = convOrderedByEnd.isEmpty()//TOREMOVE: set nextend to the end of the first element of convOrderedByEnd (or MAX_VALUE IF IT IS EMPTY)
                    ? Integer.MAX_VALUE
                    : convOrderedByEnd.get(0).getEndSite();

            if (nextStart < nextEnd) {//TOREMOVE: if the next event is a start (i.e. a conversion is added)
                nActive += 1;//TOREMOVE: increment the nb. of Active conversions

                if (nActive == 1) {//TOREMOVE: if the nb. of Active conversions is 1 (i.e. we just entered a new active region since we started from 0): we initialise the currentMergedConv
                    currentMergedConv = convOrderedByStart.get(0).getCopy();//TOREMOVE: copy the current conversion into the currentMergedConv
                    currentMergedConv.acgIndex = convOrderedByStart.get(0).acgIndex;//TOREMOVE: get acgIndex (that should be the same for all conversions)
                    mergedConvCount = 1;//TOREMOVE: initialise convcount to 1
                    mergedConvHeight1.clear();//TOREMOVE: reinitialize merged conv heights
                    mergedConvHeight2.clear();
                    mergedNodes2.clear();//TOREMOVE: reinitialize merged conv node2s
                    mergedConvHeight1.add(currentMergedConv.getHeight1());
                    mergedConvHeight2.add(currentMergedConv.getHeight2());
                    mergedNodes2.add(currentMergedConv.getNode2());
                } else {//TOREMOVE: if the nb. of Active conversions > 1 (i.e. we were already in an active region): we just increment the conv count and merged heights
                    mergedConvCount += 1;
                    mergedConvHeight1.add(convOrderedByStart.get(0).getHeight1());
                    mergedConvHeight2.add(convOrderedByStart.get(0).getHeight2());
                    mergedNodes2.add(convOrderedByStart.get(0).getNode2());
                }

                convOrderedByStart.remove(0);//TOREMOVE: remove the conversion from the convOrderedByStart list since we have included it

            } else {//TOREMOVE: if the next event is an end (i.e. a conversion is removed)
                nActive -= 1;//TOREMOVE: decrement the nb. of Active conversions

                if (nActive == 0 ) {//TOREMOVE: if the nb. of Active conversions is 0 (i.e. we just left an active region): we need to wrap up the current MergedConv and add it to the list
                    assert currentMergedConv != null;//TOREMOVE: check if the currentMergedConv is not null which shouldn't be the case
                    currentMergedConv.setEndSite(nextEnd);//TOREMOVE: set the currentMergedConv end site to the last end site we encountered before leaving the active region
                    //we get the frequency of each donor node in the conversion summary (retaining only nodes that are in the MCC CF). This is only useful for the receiverBranchMode. //TODO: receiverBranchMode
                    Map<Node, Integer> node2counts = new HashMap<>();
                    for (Node node2 : mergedNodes2) {
                        Integer count = node2counts.get(node2);
                        node2counts.put(node2, count != null ? count+1 : 1);
                    }
                    //select the most frequent node2 (or a random one among most frequent) /TODO: receiverBranchMode
                    int maxNode2count = 0;
                    Node selectedNode = null;
                    for (Node node2 : node2counts.keySet()){
                        if (node2counts.get(node2) > maxNode2count){
                            maxNode2count = node2counts.get(node2);
                            selectedNode = node2;
                        } else if (node2counts.get(node2) == maxNode2count){
                            selectedNode = Randomizer.nextBoolean() ? node2 : selectedNode;
                        }
                    }
                    //we get the sum of node heights (considering only the selected node2 in the case of height2s /TODO: receiverBranchMode
                    double sumSelectedHeights1 = 0;
                    double sumSelectedHeights2 = 0;
                    for (int i = 0; i <  mergedConvCount; i++ ){
                        if (mergedNodes2.get(i).equals(selectedNode))
                            sumSelectedHeights2 += mergedConvHeight2.get(i);
                        sumSelectedHeights1 += mergedConvHeight1.get(i);
                    }
                    //we set currentMergedConv heights and node2 and add it to the mergedList
                    currentMergedConv.setHeight1(sumSelectedHeights1 / mergedConvCount);//TOREMOVE: set the heights of the merged conversion to the mean of all included conversions by dividing the sum of heights by the conv count
                    currentMergedConv.setHeight2(sumSelectedHeights2 / maxNode2count);
                    currentMergedConv.setNode2(selectedNode);
                    mergedList.add(currentMergedConv);//TOREMOVE: add the currentMergedConv to the list
                    if (firstStep) {//TOREMOVE: circulargenome adjustment
                        numFirst = mergedConvCount;
                        firstStep = false;
                    }
                    currentMergedConv = null;//TOREMOVE: circulargenome adjustment: why necessary?
                }

                convOrderedByEnd.remove(0);//TOREMOVE: remove the conversion from the convOrderedByEnd list
            }

            if (convOrderedByEnd.isEmpty() && (nextEnd >= minOverlapStart)) {//TOREMOVE: circulargenome adjustment: if the last conversion the list was overlapping the first conversion overlapping the origin
                if (mergedList.size() > 0 && currentMergedConv != null) {//TOREMOVE: in case the list merged conversion is not empty and the current merged conversion is not null
                    mergedList.get(0).setStartSite(currentMergedConv.getStartSite());
                    //mergedList.get(0).setHeight1((mergedConvHeight1 + mergedList.get(0).getHeight1()*numFirst)/(mergedConvCount+numFirst) );
                    //mergedList.get(0).setHeight2((mergedConvHeight2 + mergedList.get(0).getHeight2()*numFirst)/(mergedConvCount+numFirst) );
                }
                break;
            }
        }

        return mergedList;
    }

    /**
     * Determine contiguous regions on specified locus where the fraction of
     * ACGs having a conversion active is greater than the given threshold.
     *
     * @param from BitSet representing source clade
     * @param to BitSet representing destination clade
     * @param locus locus to consider
     * @param threshold minimum fraction of sampled conversions included
     * @return List of regions
     */
    public List<ConversionSummary> getConversionSummaries(BitSet from, BitSet to,
                                                          Locus locus,
                                                          int nACGs,
                                                          double threshold) {

        BitSetPair bsPair = new BitSetPair(from, to);

        List<ConversionSummary> convSummaryList = new ArrayList<>();

        // Return empty list if no conversions meet the criteria.
        if (!conversionLists.containsKey(bsPair)
                || !conversionLists.get(bsPair).containsKey(locus))
            return convSummaryList;

        int thresholdCount = (int)Math.ceil(nACGs*threshold);

        List<Conversion> convOrderedByStart = new ArrayList<>();
        convOrderedByStart.addAll(conversionLists.get(bsPair).get(locus));
        convOrderedByStart.sort((Conversion o1, Conversion o2) ->
                o1.getStartSite() - o2.getStartSite());

        List<Conversion> convOrderedByEnd = new ArrayList<>();
        convOrderedByEnd.addAll(conversionLists.get(bsPair).get(locus));
        convOrderedByEnd.sort((Conversion o1, Conversion o2) ->
                o1.getEndSite() - o2.getEndSite());

        List<Conversion> activeConversions = new ArrayList<>();
        ConversionSummary conversionSummary = null;

        BitSet includedACGindices = new BitSet();

        //TODO: check adjustment (circular genome)

        int numOverlap = 0;
        for (Conversion conv : convOrderedByStart) {
            if (conv.getEndSite() < conv.getStartSite()) {
                activeConversions.add(conv);
                includedACGindices.set(conv.acgIndex);
                numOverlap += 1;
                //convOrderedByStart.remove(conv);
            }
        }

        int overlapStartBound = Integer.MAX_VALUE;
        boolean overlapRegion = false;
        if (!activeConversions.isEmpty() && activeConversions.size() >= thresholdCount) {
            overlapStartBound = activeConversions.get(thresholdCount > 0 ? thresholdCount - 1 : 0).getStartSite();
            overlapRegion = true;
            conversionSummary = new ConversionSummary();
            convSummaryList.add(conversionSummary);
            conversionSummary.addConvs(activeConversions);
        }
        int maxEndSite = convOrderedByEnd.get(convOrderedByEnd.size() -1).getEndSite();

        while (!convOrderedByStart.isEmpty() || !convOrderedByEnd.isEmpty()) {

            int nextStart = convOrderedByStart.isEmpty()
                    ? Integer.MAX_VALUE
                    : convOrderedByStart.get(0).getStartSite();

            int nextEnd = convOrderedByEnd.isEmpty()
                    ? Integer.MAX_VALUE
                    : convOrderedByEnd.get(0).getEndSite();

            if (nextStart < nextEnd) {
                activeConversions.add(convOrderedByStart.get(0));

                if (activeConversions.size() >= thresholdCount) {
                    if ( conversionSummary == null) {
                        conversionSummary = new ConversionSummary();
                        convSummaryList.add(conversionSummary);
                        conversionSummary.addConvs(activeConversions);

                        includedACGindices.clear();
                        for (Conversion conv : activeConversions)
                            includedACGindices.set(conv.acgIndex);
                    } else {
                        if (convOrderedByStart.get(0).getStartSite() <= convOrderedByStart.get(0).getEndSite()) {
                            conversionSummary.addConv(convOrderedByStart.get(0));
                            includedACGindices.set(convOrderedByStart.get(0).acgIndex);
                        }
                    }
                }
                convOrderedByStart.remove(0);
            } else {
                //TODO: check adjustment (circular genome)
                if (overlapRegion && nextEnd > overlapStartBound && nextEnd == maxEndSite) {
                    if (convSummaryList.size() > 1) {
                        convSummaryList.remove(conversionSummary);
                        convSummaryList.get(0).mergeConvSum(conversionSummary);
                        convSummaryList.get(0).nIncludedACGs += conversionSummary.nIncludedACGs;
                    }
                    conversionSummary = null;
                }
                activeConversions.remove(convOrderedByEnd.get(0));
                if (activeConversions.size() == thresholdCount-1) {
                    assert conversionSummary != null;
                    conversionSummary.nIncludedACGs = includedACGindices.cardinality();
                    conversionSummary = null;
                }
                convOrderedByEnd.remove(0);
            }
        }
        if (conversionSummary != null && thresholdCount > 0) {
            convSummaryList.remove(conversionSummary);
        }
        return convSummaryList;
    }

    /**
     * @return list of maps specifying gene flow between clades.
     */
    public List<Map<BitSet,Map<BitSet,Long>>> getGeneFlowMaps() {
        return geneFlow;
    }

    /**
     * Apply a function to each sub-clade.
     *
     * @param node MRCA of clade
     * @param function function to apply. Given sub-clade parent node
     *                 and bitset as arguments.
     * @return BitSet representing clade.
     */
    public BitSet applyToClades(Node node, BiFunction<Node, BitSet, Void> function) {
        BitSet bits = new BitSet();

        if (node.isLeaf()) {
            bits.set(2 * getTaxonIndex(node));
        } else {
            for (Node child : node.getChildren())
                bits.or(applyToClades(child, function));
        }

        function.apply(node, bits);

        return bits;
    }


        /**
         * Get bitset corresponding to a clade
         *
         * @param node MRCA of clade
         * @return BitSet representing clade.
         */
        // this method was created for the receiverBranchMode in order to get the node2 BitSets in the conversion summaries  //TODO: receiverBranchMode
        public BitSet getBitSet(Node node) {
            BitSet bitset = new BitSet();

            if (node.isLeaf()) {
                bitset.set(2 * node.getNr());
            } else {
                for (Node child : node.getChildren())
                    bitset.or(getBitSet(child));
            }

            return bitset;
        }



    /**
     * Class representing an ordered pair of BitSets.
     */
    protected class BitSetPair {
        public BitSet from, to;

        public BitSetPair(BitSet from, BitSet to) {
            this.from = from;
            this.to = to;
        }

        public BitSetPair(Conversion conv) {
            this.from = bitSets[conv.getNode1().getNr()];
            this.to = bitSets[conv.getNode2().getNr()];
        }

        @Override
        public boolean equals(Object o) {
            if (this == o) return true;
            if (o == null || getClass() != o.getClass()) return false;

            BitSetPair that = (BitSetPair) o;

            return from.equals(that.from) && to.equals(that.to);

        }

        @Override
        public int hashCode() {
            int result = from.hashCode();
            result = 31 * result + to.hashCode();
            return result;
        }

        @Override
        public String toString() {
            return from.toString() + " -> " + to.toString();
        }
    }

    /**
     * Class representing a summary of similar conversions between two
     * points in the summarized clonal frame.
     */
    public class ConversionSummary {

        List<Double> height1s = new ArrayList<>();
        List<Double> height2s = new ArrayList<>();
        List<Integer> startSites = new ArrayList<>();
        List<Integer> ends = new ArrayList<>();
        List<BitSet> node2s = new ArrayList<>();//we also store node2s in the conversion summaries (as BitSets to be consistent with the SummarizeConversion method) //TODO: receiverBranchMode

        public int nIncludedACGs = 0;

        /**
         * Add metrics associated with given conversion to summary.
         *
         * @param conv conversion
         */
        public void addConv(Conversion conv) {
            height1s.add(conv.getHeight1());
            height2s.add(conv.getHeight2());
            startSites.add(conv.getStartSite());
            ends.add(conv.getEndSite());
            node2s.add(getBitSet(conv.getNode2()));
        }

        /**
         * Add metrics associated with each of the conversions in the
         * given list to the summary.
         *
         * @param convs list of conversions
         */
        public void addConvs(List<Conversion> convs) {
            for (Conversion conv : convs)
                addConv(conv);
        }

        /**
         * Merge metrics associated in given conversion summary with this summary.
         *
         //* @param ConversionSummary convSum
         */
        public void mergeConvSum(ConversionSummary convSum) {
            height1s.addAll(convSum.height1s);
            height2s.addAll(convSum.height2s);
            startSites.addAll(convSum.startSites);
            ends.addAll(convSum.ends);
            node2s.addAll(convSum.node2s);
        }

        /**
         * @return number of conversions included in summary.
         */
        public int summarizedConvCount() {
            return height1s.size();
        }
    }
}
