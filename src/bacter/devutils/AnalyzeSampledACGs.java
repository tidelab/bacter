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
import bacter.Locus;
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

import java.io.File;
import java.io.IOException;
import java.io.PrintStream;
import java.util.*;

/**
 * @author Arthur Kocher
 */
public class AnalyzeSampledACGs {

    private static class Options {
        double burninPerc = 20.0;
        File logFile, outFile;
    }

    public static void printUsageAndExit(int exitCode) {
        System.out.println("Usage: AnalyzeSampledACGs [-burnin b] log.trees output_name");
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

        if (args.length - i < 2)
            printUsageAndExit(0);

        options.logFile = new File(args[i++]);
        options.outFile = new File(args[i++]);

        return options;
    }


    /**
     //     * Print info on each conversion in the log file
     //     */


    //main

    public static void main(String[] args) throws IOException {

        Options options = processArguments(args);

        // Set up ARG log file reader
        ACGLogReader logReader;
        logReader = new BacterACGLogReader(options.logFile, options.burninPerc);


        // Compute and write summary statistics to output file

        try (PrintStream ps = new PrintStream(options.outFile)) {
            ps.println("iter" + "\t" + "locus" + "\t" + "conv.id" + "\t" + "start.site" + "\t" + "end.site");
            int i = 0;
            for (ConversionGraph acg : logReader) {

                for (Locus locus : acg.getConvertibleLoci()) {

                    int j = 0;

                    for (Conversion conv : acg.getConversions(locus)) {

                        ps.println(i + "\t" + locus.getID() + "\t" + j + "\t" +
                                conv.getStartSite() + "\t" + conv.getEndSite());

                        j++;

                    }

                }

                i++;

            }
        }
    }
}

