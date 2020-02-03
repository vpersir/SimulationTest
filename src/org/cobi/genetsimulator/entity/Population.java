/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.cobi.genetsimulator.entity;

import cern.colt.bitvector.BitVector;

import cern.colt.list.IntArrayList;
import cern.colt.matrix.DoubleMatrix2D;
import cern.colt.matrix.impl.DenseDoubleMatrix2D;
import cern.jet.random.sampling.RandomSampler;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.LinkedList;
import java.util.List;
import org.apache.commons.math.random.MersenneTwister;
import org.apache.commons.math.random.RandomData;
import org.apache.commons.math.random.RandomDataImpl;
import org.cobi.genetsimulator.controller.LiabilityThresholdModel.PredisposingSNP;
import org.cobi.genetsimulator.controller.MultiLocusLiabilityThresholdModel;
import org.cobi.genetsimulator.controller.MultiLocusLiabilityThresholdModel.LiabilityDiseaseSNP;

import org.cobi.genetsimulator.controller.MultilocusLogitModel;
import org.cobi.genetsimulator.controller.MultilocusLogitModel.LogitDiseaseSNP;
import org.cobi.genetsimulator.controller.PopuStatSummarizer;
import org.cobi.genetsimulator.controller.RischMultipLociModel;
import org.cobi.util.stat.MultiNormalRandGenerator;
import umontreal.ssj.randvar.NormalGen;
import umontreal.ssj.randvarmulti.MultinormalCholeskyGen;
import umontreal.ssj.rng.MT19937;
import umontreal.ssj.rng.WELL607;
 
/**
 *
 * @author mxli
 */
public class Population {

    List<Individual> allIndiv;
    MultiNormalRandGenerator multiNormGenerator;
    RandomData randGenerator;

    public Population() {
        this.allIndiv = Collections.synchronizedList(new LinkedList<Individual>());
        //this.allIndiv =Collections.synchronizedList(new ArrayList<Individual>());
        //this.allIndiv = new ArrayList<Individual>();
        this.multiNormGenerator = new MultiNormalRandGenerator();
        this.randGenerator = new RandomDataImpl(new MersenneTwister());
    }

    public synchronized List<Individual> getAllIndiv() {
        return allIndiv;
    }

    public void setAllIndiv(List<Individual> allIndiv) {
        this.allIndiv = allIndiv;
    }

    public Population(int indivSize) {
        allIndiv = new ArrayList<Individual>(indivSize);

        for (int i = 0; i < indivSize; i++) {
            Individual indi = new Individual();
            indi.setFamilyID(String.valueOf(i));
            indi.setIndividualID("0");
            indi.setDadID("0");
            indi.setMomID("0");
            indi.setGender(1);
            indi.markerGtySet = new StatusGtySet();
            indi.traitGtySet = new StatusGtySet();
            allIndiv.add(indi);
           // indi.getTraits().add(null);
        }
        this.multiNormGenerator = new MultiNormalRandGenerator();
        this.randGenerator = new RandomDataImpl(new MersenneTwister());
    }

    public void readReferenceHapMapPhase() {
    }

    public boolean calculateGenotypeCovarianceMatrix(DoubleMatrix2D jointProb) throws Exception {
        multiNormGenerator.setJointProb(jointProb);
        return multiNormGenerator.exploreCovarMatrix();
    }

    public void allocateMarkerGenotypeSpaces(int lociNum) throws Exception {
        int lociNumLess = lociNum - 1;
        int size = allIndiv.size();
        for (int i = 0; i < size; i++) {
            Individual indi = allIndiv.get(i);
            indi.markerGtySet.existence = new BitVector(lociNum);
            indi.markerGtySet.existence.replaceFromToWith(0, lociNumLess, true);
            indi.markerGtySet.paternalChrom = new BitVector(lociNum);
            indi.markerGtySet.paternalChrom.replaceFromToWith(0, lociNumLess, false);
            indi.markerGtySet.maternalChrom = new BitVector(lociNum);
            indi.markerGtySet.maternalChrom.replaceFromToWith(0, lociNumLess, false);
        }
    }

    public void allocateMarkerTraitGenotypeSpaces(int markerLociNum, int traitLociNum) throws Exception {
        int markerLociNumLess = markerLociNum - 1;
        int traitLociNumLess = traitLociNum - 1;
        int size = allIndiv.size();
        if (markerLociNum > 0) {
            for (int i = 0; i < size; i++) {
                Individual indi = allIndiv.get(i);
                indi.markerGtySet.existence = new BitVector(markerLociNum);
                indi.markerGtySet.existence.replaceFromToWith(0, markerLociNumLess, true);
                indi.markerGtySet.paternalChrom = new BitVector(markerLociNum);
                indi.markerGtySet.paternalChrom.replaceFromToWith(0, markerLociNumLess, false);
                indi.markerGtySet.maternalChrom = new BitVector(markerLociNum);
                indi.markerGtySet.maternalChrom.replaceFromToWith(0, markerLociNumLess, false);
            }
        }

        if (traitLociNum > 0) {
            for (int i = 0; i < size; i++) {
                Individual indi = allIndiv.get(i);
                indi.traitGtySet.existence = new BitVector(traitLociNum);
                indi.traitGtySet.existence.replaceFromToWith(0, traitLociNumLess, true);
                indi.traitGtySet.paternalChrom = new BitVector(traitLociNum);
                indi.traitGtySet.paternalChrom.replaceFromToWith(0, traitLociNumLess, false);
                indi.traitGtySet.maternalChrom = new BitVector(traitLociNum);
                indi.traitGtySet.maternalChrom.replaceFromToWith(0, traitLociNumLess, false);
            }
        }

    }

    /*
    //strange I found CorrelatedRandomVectorGenerator does not work for a matrix with over 6 variables
    public void simulateDependentGenotypes1() throws Exception {
    DoubleMatrix2D covarianceMatrix = multiNormGenerator.getCovarianceMatrix();
    double[] quantileProb = multiNormGenerator.getQuantileProb();
    int lociNum = quantileProb.length;
    int lociNumLess = lociNum - 1;
    final double PRECISION = 1.0e-8;
    double[] mean = new double[lociNum];
    Arrays.fill(mean, 0.0);
    int indivSize = allIndiv.size();
    
    CorrelatedRandomVectorGenerator sg = new CorrelatedRandomVectorGenerator(mean, covarianceMatrix, PRECISION, new GaussianRandomGenerator(new MersenneTwister()));
    NormalGen ng = new NormalGen(new MT19937(new WELL607()));
    
    // System.out.println(mean.toString());
    // System.out.println(covarianceMatrix.toString());
    
    for (int i = 0; i < indivSize; i++) {
    double[] samples = sg.nextVector();
    allIndiv.get(i).markerGtySet.paternalChrom.replaceFromToWith(0, lociNumLess, false);
    for (int j = 0; j < lociNum; j++) {
    if (samples[j] <= quantileProb[j]) {
    allIndiv.get(i).markerGtySet.paternalChrom.putQuick(j, true);
    }
    }
    samples = sg.nextVector();
    System.out.println(samples[lociNum - 1]);
    allIndiv.get(i).markerGtySet.maternalChrom.replaceFromToWith(0, lociNumLess, false);
    for (int j = 0; j < lociNum; j++) {
    if (samples[j] <= quantileProb[j]) {
    allIndiv.get(i).markerGtySet.maternalChrom.putQuick(j, true);
    }
    }
    }
    }
     */
    public void simulateDependentGenotypes() throws Exception {
        DoubleMatrix2D covarianceMatrix = multiNormGenerator.getCovarianceMatrix();
        double[] quantileProb = multiNormGenerator.getQuantileProb();
        int lociNum = quantileProb.length;
        int lociNumLess = lociNum - 1;
        final double PRECISION = 1.0e-8;
        double[] mean = new double[lociNum];
        Arrays.fill(mean, 0.0);
        int indivSize = allIndiv.size();

        NormalGen ng = new NormalGen(new MT19937(new WELL607()));

        // System.out.println(mean.toString());
        //  System.out.println(covarianceMatrix.toString());
        MultinormalCholeskyGen sg = new MultinormalCholeskyGen(ng, mean, covarianceMatrix); //How about the alternative one MultinormalPCAGen?
        double[] samples = new double[lociNum];
        for (int i = 0; i < indivSize; i++) {
            sg.nextPoint(samples);
            allIndiv.get(i).markerGtySet.paternalChrom.replaceFromToWith(0, lociNumLess, false);
            for (int j = 0; j < lociNum; j++) {
                if (samples[j] <= quantileProb[j]) {
                    allIndiv.get(i).markerGtySet.paternalChrom.putQuick(j, true);
                }
            }
            sg.nextPoint(samples);
            // System.out.println(samples[lociNum - 1]);
            allIndiv.get(i).markerGtySet.maternalChrom.replaceFromToWith(0, lociNumLess, false);
            for (int j = 0; j < lociNum; j++) {
                if (samples[j] <= quantileProb[j]) {
                    allIndiv.get(i).markerGtySet.maternalChrom.putQuick(j, true);
                }
            }
        }
    }

    public void simulateDependentGenotypesByBlock(int blockNum) throws Exception {
        DoubleMatrix2D covarianceMatrix1 = multiNormGenerator.getCovarianceMatrix();
        double[] quantileProb = multiNormGenerator.getQuantileProb();
        int lociNum = quantileProb.length;

        final double PRECISION = 1.0e-8;
        double[] mean = new double[lociNum];
        Arrays.fill(mean, 0.0);
        int indivSize = allIndiv.size();
        int[] seeds = new int[19];
        for (int i = 0; i < seeds.length; i++) {
            seeds[i] = (int) (Math.random() * 1000000);
           // System.out.println(seeds[i]);
        }
        WELL607 we = new WELL607();
        we.setSeed(seeds);
 
        NormalGen ng = new NormalGen(new MT19937(we));       
        DenseDoubleMatrix2D covarianceMatrix = new DenseDoubleMatrix2D(lociNum, lociNum);
        for (int i = 0; i < lociNum; i++) {
            for (int j = 0; j < lociNum; j++) {
                covarianceMatrix.setQuick(i, j, covarianceMatrix1.getQuick(i, j));
            }
        }
        // System.out.println(mean.toString());
        // System.out.println(covarianceMatrix.toString());
        MultinormalCholeskyGen sg = new MultinormalCholeskyGen(ng, mean, covarianceMatrix); //How about the alternative one MultinormalPCAGen?
        double[] samples = new double[lociNum];
        //CorrelatedRandomVectorGenerator sg = new CorrelatedRandomVectorGenerator(mean, covarianceMatrix, PRECISION, new GaussianRandomGenerator(new MersenneTwister()));
        for (int b = 0; b < blockNum; b++) {
            for (int i = 0; i < indivSize; i++) {
                sg.nextPoint(samples);
                allIndiv.get(i).markerGtySet.paternalChrom.replaceFromToWith(b * lociNum, (b + 1) * lociNum - 1, false);
                for (int j = 0; j < lociNum; j++) {
                    if (samples[j] <= quantileProb[j]) {
                        allIndiv.get(i).markerGtySet.paternalChrom.putQuick(b * lociNum + j, true);
                    }
                }
                sg.nextPoint(samples);
                allIndiv.get(i).markerGtySet.maternalChrom.replaceFromToWith(b * lociNum, (b + 1) * lociNum - 1, false);
                for (int j = 0; j < lociNum; j++) {
                    if (samples[j] <= quantileProb[j]) {
                        allIndiv.get(i).markerGtySet.maternalChrom.putQuick(b * lociNum + j, true);
                    }
                }
            }
        }
    }

    public void simulateIndependentMarkerGenotypes(double[] markerAlleleFreqs) throws Exception {
        int lociNum = markerAlleleFreqs.length;
        int lociNumLess = lociNum - 1;
        int indivSize = allIndiv.size();
        double prob;

        for (int i = 0; i < indivSize; i++) {
            allIndiv.get(i).markerGtySet.paternalChrom.replaceFromToWith(0, lociNumLess, false);
            for (int j = 0; j < lociNum; j++) {
                prob = randGenerator.nextUniform(0, 1);
                if (prob <= markerAlleleFreqs[j]) {
                    allIndiv.get(i).markerGtySet.paternalChrom.putQuick(j, true);
                }
            }

            allIndiv.get(i).markerGtySet.maternalChrom.replaceFromToWith(0, lociNumLess, false);
            for (int j = 0; j < lociNum; j++) {
                prob = randGenerator.nextUniform(0, 1);
                if (prob <= markerAlleleFreqs[j]) {
                    allIndiv.get(i).markerGtySet.maternalChrom.putQuick(j, true);
                }
            }
        }
    }

    public void simulatePhenotypeMultiLiabilityThresholdModel(MultiLocusLiabilityThresholdModel multiLiabModel) throws Exception {
        List<LiabilityDiseaseSNP> liabiSNPs = multiLiabModel.getLiabiSNPs();

        int lociNum = liabiSNPs.size();
        int lociNumLess = lociNum - 1;
        int indivSize = allIndiv.size();
        double prob;
        double accumulatedLiabilities = 0.0;
        double residualVariance = Math.sqrt(1 - multiLiabModel.getAccumulatedVariance());

        double overallThreshold = multiLiabModel.getOverallThreshold();
        int[] suscepLoci = new int[lociNum];
        for (int j = 0; j < lociNum; j++) {
            suscepLoci[j] = liabiSNPs.get(j).getLDMarkerPosition();
        }

        for (int i = 0; i < indivSize; i++) {
            allIndiv.get(i).traitGtySet.paternalChrom.replaceFromToWith(0, lociNumLess, false);
            for (int j = 0; j < lociNum; j++) {
                prob = randGenerator.nextUniform(0, 1);
                //assume M and A alleles are encoded as "true".
                if (allIndiv.get(i).markerGtySet.paternalChrom.getQuick(suscepLoci[j])) {
                    if (prob <= liabiSNPs.get(j).probAcM) {
                        allIndiv.get(i).traitGtySet.paternalChrom.putQuick(j, true);
                    }
                } else if (prob <= liabiSNPs.get(j).probAcm) {
                    allIndiv.get(i).traitGtySet.paternalChrom.putQuick(j, true);
                }
            }

            allIndiv.get(i).traitGtySet.maternalChrom.replaceFromToWith(0, lociNumLess, false);
            for (int j = 0; j < lociNum; j++) {
                prob = randGenerator.nextUniform(0, 1);
                //assume M and A alleles are encoded as "true".
                if (allIndiv.get(i).markerGtySet.maternalChrom.getQuick(suscepLoci[j])) {
                    if (prob <= liabiSNPs.get(j).probAcM) {
                        allIndiv.get(i).traitGtySet.maternalChrom.putQuick(j, true);
                    }
                } else if (prob <= liabiSNPs.get(j).probAcm) {
                    allIndiv.get(i).traitGtySet.maternalChrom.putQuick(j, true);
                }
            }
            accumulatedLiabilities = 0.0;

            for (int j = 0; j < lociNum; j++) {
                if (allIndiv.get(i).traitGtySet.paternalChrom.getQuick(j) && allIndiv.get(i).traitGtySet.maternalChrom.getQuick(j)) {
                    accumulatedLiabilities += liabiSNPs.get(j).liabilityAA;
                } else if (!allIndiv.get(i).traitGtySet.paternalChrom.getQuick(j) && !allIndiv.get(i).traitGtySet.maternalChrom.getQuick(j)) {
                    accumulatedLiabilities += liabiSNPs.get(j).liabilityaa;
                } else {
                    accumulatedLiabilities += liabiSNPs.get(j).liabilityAa;
                }
            }

            prob = randGenerator.nextGaussian(accumulatedLiabilities, residualVariance);
            if (prob >= overallThreshold) {
                allIndiv.get(i).setAffectedStatus(2);
            } else {
                allIndiv.get(i).setAffectedStatus(1);
            }
        }
    }

    public void simulatePhenotypRischMultipLociModel(RischMultipLociModel rischModel) throws Exception {
        int lociNum = rischModel.getLociNum();

        int indivSize = allIndiv.size();
        double prob;
        List<DiseaseSNP> liabiSNPs = rischModel.getDiseaseSNPs();
        int[] suscepLociID = new int[lociNum];
        for (int j = 0; j < lociNum; j++) {
            suscepLociID[j] = liabiSNPs.get(j).getLDMarkerPosition();
        }
        boolean[] suscepLociLabel = new boolean[lociNum];
        for (int j = 0; j < lociNum; j++) {
            suscepLociLabel[j] = liabiSNPs.get(j).riskAlleleLable;
        }
        double[] jointPenetrance = rischModel.getJointGenotypePenerances();
        int riskAlleleNum = 0;
        for (int i = 0; i < indivSize; i++) {
            prob = randGenerator.nextUniform(0, 1);
            riskAlleleNum = 0;

            for (int j = 0; j < lociNum; j++) {
                //directly use the markerGtySet to decide the disease probaility
                if (allIndiv.get(i).markerGtySet.paternalChrom.getQuick(suscepLociID[j]) == suscepLociLabel[j]) {
                    riskAlleleNum++;
                }
                if (allIndiv.get(i).markerGtySet.maternalChrom.getQuick(suscepLociID[j]) == suscepLociLabel[j]) {
                    riskAlleleNum++;
                }
            }
            if (prob <= jointPenetrance[riskAlleleNum]) {
                allIndiv.get(i).setAffectedStatus(2);
            } else {
                allIndiv.get(i).setAffectedStatus(1);
            }
        }
    }

    public void simulatePhenotypeMultilocusLogitModel(MultilocusLogitModel multiLogitModel) throws Exception {
        List<LogitDiseaseSNP> snpList = multiLogitModel.getSnpList();
        double alpha = multiLogitModel.getAlpha();
        int suscepNum = snpList.size();
        int[] suscepLoci = new int[suscepNum];
        boolean[] riskAlleles = new boolean[suscepNum];
        int indivSize = allIndiv.size();
        for (int j = 0; j < suscepNum; j++) {
            suscepLoci[j] = snpList.get(j).getLDMarkerPosition();
            riskAlleles[j] = snpList.get(j).isRiskAlleleLable();
        }
        double alphaBeta;
        double jointGtyPenerance;
        for (int i = 0; i < indivSize; i++) {
            double prob = randGenerator.nextUniform(0, 1);
            StatusGtySet gty = allIndiv.get(i).markerGtySet;
            alphaBeta = alpha;
            for (int j = 0; j < suscepNum; j++) {
                if (gty.paternalChrom.getQuick(suscepLoci[j]) == riskAlleles[j] && gty.maternalChrom.getQuick(suscepLoci[j]) == riskAlleles[j]) {
                    alphaBeta += snpList.get(j).getBeta2();
                } else if (gty.paternalChrom.getQuick(suscepLoci[j]) != riskAlleles[j] && gty.maternalChrom.getQuick(suscepLoci[j]) != riskAlleles[j]) {
                    //not change alphaBeta
                } else {
                    alphaBeta += snpList.get(j).getBeta1();
                }
            }
            alphaBeta = Math.exp(alphaBeta);
            jointGtyPenerance = 1 - 1 / (1 + alphaBeta);

            if (prob < jointGtyPenerance) {
                allIndiv.get(i).setAffectedStatus(2);
            } else {
                allIndiv.get(i).setAffectedStatus(1);
            }
        }

        PopuStatSummarizer popuSum = new PopuStatSummarizer();
        popuSum.summarizePhenotypePropertyforLogitDiseaseSNP(allIndiv, snpList);
    }

    public void simulatePhenotypeRandom(double prevlance) throws Exception {
        int indivSize = allIndiv.size();
        int caseNum = (int) (indivSize * prevlance);
        IntArrayList caseIndexes = new IntArrayList();
        for (int i = 0; i < indivSize; i++) {
            allIndiv.get(i).setAffectedStatus(1);
            caseIndexes.add(i);
        }
        //http://acs.lbl.gov/software/colt/api/cern/jet/random/sampling/RandomSampler.html#RandomSampler(long, long, long, cern.jet.random.engine.RandomEngine)       
        long[] caseindexes = new long[caseNum];
        RandomSampler.sample(caseNum, caseIndexes.size(), caseNum, 0, caseindexes, 0, new cern.jet.random.engine.MersenneTwister(new java.util.Date()));
        for (int i = 0; i < caseNum; i++) {
            allIndiv.get(caseIndexes.getQuick((int) caseindexes[i])).setAffectedStatus(2);
        }


    }

    public void simulatePhenotypeLiabilityThresholdModel(List<PredisposingSNP> snpList) throws Exception {
        int suscepNum = snpList.size();
        int[] suscepLoci = new int[suscepNum];
        boolean[] riskAlleles = new boolean[suscepNum];
        int indivSize = allIndiv.size();
        for (int j = 0; j < suscepNum; j++) {
            suscepLoci[j] = snpList.get(j).getLDMarkerPosition();
            riskAlleles[j] = snpList.get(j).isRiskAlleleLable();
        }

        StringBuffer susceptGtyLabel = new StringBuffer();
        susceptGtyLabel.setLength(suscepNum);

        for (int i = 0; i < indivSize; i++) {
            double prob = randGenerator.nextUniform(0, 1);
            StatusGtySet gty = allIndiv.get(i).markerGtySet;
            for (int j = 0; j < suscepNum; j++) {
                if (gty.paternalChrom.getQuick(suscepLoci[j]) == riskAlleles[j] && gty.maternalChrom.getQuick(suscepLoci[j]) == riskAlleles[j]) {
                    susceptGtyLabel.setCharAt(j, '2');
                } else if (gty.paternalChrom.getQuick(suscepLoci[j]) != riskAlleles[j] && gty.maternalChrom.getQuick(suscepLoci[j]) != riskAlleles[j]) {
                    susceptGtyLabel.setCharAt(j, '0');
                } else {
                    susceptGtyLabel.setCharAt(j, '1');
                }
            }
            /*
            if (prob < jointGtyPenerance[LocalNumber.xSystem2decimal(susceptGtyLabel, 3)]) {
            allIndiv.get(i).setAffectedStatus(2);
            } else {
            allIndiv.get(i).setAffectedStatus(1);
            }*/
        }

        PopuStatSummarizer popuSum = new PopuStatSummarizer();
        //  popuSum.summarizePhenotypeProperty(allIndiv, snpList);
    }

    public void simulatePhenotypeContinuousTraitPolygenicAdditiveModel(int[] suscepLoci, boolean[] riskAlleles, double alaph, double ve) throws Exception {
        int suscepNum = suscepLoci.length;
        int indivSize = allIndiv.size();
        double gtyScore = 0;
        for (int i = 0; i < indivSize; i++) {
            double prob = randGenerator.nextUniform(0, 1);
            StatusGtySet gty = allIndiv.get(i).markerGtySet;
            gtyScore = 0;
            for (int j = 0; j < suscepNum; j++) {
                if (gty.paternalChrom.getQuick(suscepLoci[j]) == riskAlleles[j] && gty.maternalChrom.getQuick(suscepLoci[j]) == riskAlleles[j]) {
                    gtyScore += 2;
                } else if (gty.paternalChrom.getQuick(suscepLoci[j]) != riskAlleles[j] && gty.maternalChrom.getQuick(suscepLoci[j]) != riskAlleles[j]) {
                } else {
                    gtyScore += 1;
                }
            }
            gtyScore = gtyScore * alaph + randGenerator.nextGaussian(0, Math.sqrt(ve));
            allIndiv.get(i).setMainTrait(new double[]{gtyScore});
        }
    }

    public void simulateTwoPhenotypeContinuousTraitPolygenicAdditiveModel(int[] suscepLoci, boolean[] riskAlleles, double alaph1, double ve1, double alaph2, double ve2) throws Exception {
        int suscepNum = suscepLoci.length;
        int indivSize = allIndiv.size();
        double gtyScore1 = 0;
        double gtyScore2 = 0;
        for (int i = 0; i < indivSize; i++) {
            double prob = randGenerator.nextUniform(0, 1);
            StatusGtySet gty = allIndiv.get(i).markerGtySet;
            gtyScore1 = 0;
            for (int j = 0; j < suscepNum; j++) {
                if (gty.paternalChrom.getQuick(suscepLoci[j]) == riskAlleles[j] && gty.maternalChrom.getQuick(suscepLoci[j]) == riskAlleles[j]) {
                    gtyScore1 += 2;
                } else if (gty.paternalChrom.getQuick(suscepLoci[j]) != riskAlleles[j] && gty.maternalChrom.getQuick(suscepLoci[j]) != riskAlleles[j]) {
                } else {
                    gtyScore1 += 1;
                }
            }
            gtyScore2 = gtyScore1;
            gtyScore1 = gtyScore1 * alaph1 + randGenerator.nextGaussian(0, Math.sqrt(ve1));
            gtyScore2 = gtyScore2 * alaph2 + randGenerator.nextGaussian(0, Math.sqrt(ve2));
            allIndiv.get(i).setMainTrait(new double[]{gtyScore1, gtyScore2});
        }
    }

    public void nullifyPoulation() {
        int size = allIndiv.size();
        for (int i = 0; i < size; i++) {
            allIndiv.set(i, null);
        }
    }
}
