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
import bacter.acgannotator.ACGAnnotator;
import bacter.model.ACGLikelihood;
import bacter.util.ACGLogReader;
import bacter.util.BacterACGLogReader;
import bacter.util.COACGLogFileReader;
import beast.evolution.alignment.Alignment;
import beast.evolution.sitemodel.SiteModel;
import beast.evolution.substitutionmodel.Frequencies;
import beast.evolution.substitutionmodel.HKY;
import beast.evolution.tree.Node;
import beast.evolution.tree.TreeUtils;
import beast.util.NexusParser;

import javax.xml.stream.XMLStreamException;
import java.io.File;
import java.io.IOException;
import java.io.PrintStream;
import java.util.*;

/**
 * @author Tim Vaughan <tgvaughan@gmail.com>
 */
public class SummarizedVsTrueACG {

    private static class Options {
        double burninPerc = 20.0;
        File logFile, truthFile, aliFile, summarizedAcgOut, outFileSummarizedClade, outFileSummarizedConv, outFileTrueConv;
    }

    public static void printUsageAndExit(int exitCode) {
        System.out.println("Usage: SummarizedVsTrueACG [-burnin b] truth.tree log.trees true_alignment_file output_summarized_ACG output_compare_summarized_clades output_compare_summarized_conversions output_compare_true_conversions");
        System.exit(exitCode);
    }

    /**
     * Process command line arguments.
     *
     * @param args list of arguments given in the command line
     * @return an Option variable
     */
    public static Options processArguments(String[] args) {

        Options options = new Options();

        int i = 0;
        while (i < args.length && args[i].startsWith("-")) {
            switch (args[i].substring(1)) {

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

    /**
     //     * Get and print info on summarized ACG clades (if they are true and their support)
     //     *
     //     * @param trueClades clades in true acg
     //     * @param summarizedACG summary of sampled ACGs
     //     * @ps print stream
     //     */
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

    /**
     //     * Get and print info on conversions in summarized ACG (if they are true, their support and the error on the converted region edges)
     //     *
     //     * @param trueACG true acg
     //     * @param trueClades clades in true acg
     //     * @param summarizedACG summary of sampled ACGs
     //     * @ps print stream
     //     */
    public static void compareSumACGConv(ConversionGraph trueACG, Clade[] trueClades,
                                         ConversionGraph summarizedACG, PrintStream ps) {

        Clade[] clades = new Clade[summarizedACG.getNodeCount()];
        getClades(clades, summarizedACG.getRoot());
        int convNum = 1;
        ps.println("ConvNb" + "\t" + "isTrue" + "\t" + "support" + "\t" + "edgeError");

        for (Conversion conv : summarizedACG.getConversions(summarizedACG.getConvertibleLoci().get(0))) {

            int isTrue = 0;
            double edgeError = .0;

            //get parent clades of the conversion
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
            double posterior = Double.parseDouble(conv.newickMetaDataMiddle.replaceAll(".*posterior=([^,]*).*", "$1"));

            for (Conversion trueConv : trueACG.getConversions(trueACG.getConvertibleLoci().get(0))) {
                Clade trueFromClade = trueClades[trueConv.getNode1().getNr()];
                Clade trueToClade = trueClades[trueConv.getNode2().getNr()];
                int trueStartSite = trueConv.getStartSite();
                int trueEndSite = trueConv.getEndSite();
                //we define a conversion as true if there is a conversion in the true ACG that have the same parent clades and an overlap of at least 1 bp in the converted region
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

    /**
     //     * Get and print info on true conversions (posterior support, error on the estimated converted region edges, length of converted region, start height, end height, likelihood gain)
     //     *
     //     * @param trueACG true acg
     //     * @param trueClades clades in true acg
     //     * @param summarizedACG summary of sampled ACGs
     //     * @param trueAlignment alignment from which ACGs where estimated
     //     * @ps print stream
     //     */

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
        ps.println("TrueConvNb" + "\t" + "supportInSummarizedACG" + "\t" + "edgeError" + "\t" + "tractLength" + "\t" + "startHeight" + "\t" + "endHeigth" + "\t" + "likelihoodGain" + "\t" + "distanceBetweenConv");

        for (Conversion trueConv : trueACG.getConversions(trueACG.getConvertibleLoci().get(0))) {
            Clade trueFromClade = trueClades[trueConv.getNode1().getNr()];
            Clade trueToClade = trueClades[trueConv.getNode2().getNr()];
            int trueStartSite = trueConv.getStartSite();
            int trueEndSite = trueConv.getEndSite();
            double edgeError = 1;
            double posterior = 0.0;
            int tractLength = trueConv.getSiteCount();
            double startHeight = trueConv.getHeight2();
            double endHeight = trueConv.getHeight1();
            //get common ancestor of parental clades
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
                convDist = trueConv.getHeight2() - trueConv.getHeight1();
            } else{
                convDist = Math.abs(commonAncestor.getHeight()-trueConv.getHeight1()) + Math.abs(commonAncestor.getHeight()-trueConv.getHeight2());
            }
            //get the gain in tree likelihood associated by this conversion (ratio of tree likelihood with our without the conversion)
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
                //we define a conversion as true if there is a conversion in the true ACG that have the same parent clades and an overlap of at least 1 bp in the converted region
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
            ps.println(convNum++ + "\t" + posterior + "\t" + edgeError + "\t" + tractLength + "\t" + startHeight + "\t" + endHeight + "\t" + likelihoodGain + "\t" + convDist);
        }
    }

    //main

    public static void main(String[] args) throws IOException {

        Options options = processArguments(args);

        // Load true ARG

        BacterACGLogReader truthReader = new BacterACGLogReader(options.truthFile, 0);
        if (truthReader.getACGCount() != 1) {
            System.out.println("Expected exactly 1 ACG in truth file. Found " +
                    truthReader.getACGCount());
            System.exit(1);
        }

        ConversionGraph trueACG = null;
        for (ConversionGraph acg : truthReader) trueACG = acg;

        // Determine clades present in truth
        assert trueACG != null;
        Clade[] trueClades;
        trueClades = new Clade[trueACG.getNodeCount()];
        getClades(trueClades, trueACG.getRoot());

        // Load alignments
        Alignment trueAlignment;
        NexusParser nexusParser = new NexusParser();
        nexusParser.parseFile(options.aliFile);
        trueAlignment = nexusParser.m_alignment;

        // Compute summary ACG
        ACGAnnotator acgAnnotator = new ACGAnnotator(options.logFile, options.summarizedAcgOut, options.burninPerc, 1);
        ConversionGraph summarizedACG = acgAnnotator.getSummarizedACG();

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
