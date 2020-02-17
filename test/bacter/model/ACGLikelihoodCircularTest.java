package bacter.model;

import bacter.Conversion;
import bacter.ConversionGraph;
import bacter.Locus;
import bacter.TestBase;
import beast.core.parameter.RealParameter;
import beast.evolution.sitemodel.SiteModel;
import beast.evolution.substitutionmodel.JukesCantor;
import beast.evolution.tree.Node;
import beast.evolution.tree.coalescent.ConstantPopulation;
import beast.util.ClusterTree;
import org.junit.Test;

import static org.junit.Assert.assertTrue;

public class ACGLikelihoodCircularTest extends TestBase {

    public ACGLikelihoodCircularTest() { }

    @Test
    public void testLikelihoodFixedData() throws Exception {

        Locus locus = new Locus("locus", getSimpleAlignment());

        // ConversionGraph
        ConversionGraph acg = new ConversionGraph();
        ClusterTree tree = new ClusterTree();
        tree.initByName(
                "clusterType", "upgma",
                "taxa", locus.getAlignment());

        acg.assignFrom(tree);
        acg.initByName("locus", locus,
                                "circularGenome", "true");

        // Site model:
        JukesCantor jc = new JukesCantor();
        jc.initByName();
        SiteModel siteModel = new SiteModel();
        siteModel.initByName(
                "substModel", jc);

        // Likelihood

        ACGLikelihood argLikelihood = new ACGLikelihood();
        argLikelihood.initByName(
                "locus", locus,
                "tree", acg,
                "siteModel", siteModel);

        acg.setEverythingDirty(true);

        //Add a single recombination event
        acg.getNodesAsArray();
        Node node1 = acg.getExternalNodes().get(0); //acg.getExternalNodes().get(0);
        Node node2 = node1.getParent();
        double height1 = 0.5*(node1.getHeight() + node1.getParent().getHeight());
        double height2 = 0.5*(node2.getHeight() + node2.getParent().getHeight());
        int startLocus = 1; //0;
        int endLocus = 3; //6;

        Conversion recomb1 = new Conversion(node1, height1, node2, height2,
                startLocus, endLocus, acg, locus);
        acg.addConversion(recomb1);

        double logP = argLikelihood.calculateLogP();

        //Define recombination on other half of circular genome
        acg.deleteConversion(recomb1);

        startLocus = 4; //7;
        endLocus = 0; //13;
        Conversion recomb2 = new Conversion(node1, height1, node2, height2,
                startLocus, endLocus, acg, locus);
        acg.addConversion(recomb2);
        //acg.initAndValidate();
        //argLikelihood.initAndValidate();

        double logPOtherHalf = argLikelihood.calculateLogP();

        double relativeDiff = Math.abs(2.0*(logPOtherHalf-logP)/(logPOtherHalf+logP));

        assertTrue(relativeDiff<1e-14);


    }

}
