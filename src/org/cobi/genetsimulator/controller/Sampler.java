/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.cobi.genetsimulator.controller;

import cern.colt.list.DoubleArrayList;
import cern.colt.list.IntArrayList;
import cern.colt.matrix.DoubleMatrix2D;
import cern.jet.random.sampling.RandomSampler;
import cern.jet.stat.Descriptive;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Iterator;
import java.util.List;
import org.apache.commons.math.random.MersenneTwister;
import org.apache.commons.math.random.RandomData;
import org.apache.commons.math.random.RandomDataImpl;
import org.cobi.genetsimulator.entity.Individual;
import org.cobi.genetsimulator.entity.IndividualComparator;
import org.cobi.genetsimulator.entity.StatusGtySet;
import umontreal.ssj.randvar.NormalGen;
import umontreal.ssj.randvarmulti.MultinormalCholeskyGen;
import umontreal.ssj.randvarmulti.MultinormalGen;
import umontreal.ssj.rng.MT19937;
import umontreal.ssj.rng.WELL607;
 

/**
 *
 * @author mxli
 */
public class Sampler {

    RandomData randGenerator;
 
    public Sampler() {
        randGenerator = new RandomDataImpl(new MersenneTwister());
    }

    public List<Individual> sampleWithoutReplacement(List<Individual> allIndiv, int caseNum, int controlNum) {
        Collections.shuffle(allIndiv);
        List<Individual> sampledIndiv = new ArrayList<Individual>();
        Iterator<Individual> iter = allIndiv.iterator();
        while (iter.hasNext()) {
            Individual indiv = iter.next();
            if (caseNum > 0 && indiv.getAffectedStatus() == 2) {
                sampledIndiv.add(indiv);
                caseNum--;
            } else if (controlNum > 0 && indiv.getAffectedStatus() == 1) {
                sampledIndiv.add(indiv);
                controlNum--;
            }
            if (caseNum == 0 && controlNum == 0) {
                break;
            }
        }
        if (caseNum != 0) {
            System.out.println("Warning! The population has less case than " + caseNum);
        }

        if (controlNum != 0) {
            System.out.println("Warning! The population has less control than " + controlNum);
        }
        return sampledIndiv;
    }

    public List<Individual> sampleWithoutReplacement(List<Individual> allIndiv, int caseNum, int controlNum, IntArrayList caseIndexes, IntArrayList controlIndexes) {
        //very time consumming for large list
        //Collections.shuffle(allIndiv);
        if (caseNum >= caseIndexes.size()) {
            System.out.println("Warning! The population has less case than " + caseNum);
        }

        if (controlNum >= controlIndexes.size()) {
            System.out.println("Warning! The population has less control than " + controlNum);
        }

        //http://acs.lbl.gov/software/colt/api/cern/jet/random/sampling/RandomSampler.html#RandomSampler(long, long, long, cern.jet.random.engine.RandomEngine)       
        long[] caseindexes = new long[caseNum];
        RandomSampler.sample(caseNum, caseIndexes.size(), caseNum, 0, caseindexes, 0, new cern.jet.random.engine.MersenneTwister(new java.util.Date()));
        List<Individual> sampledIndiv = new ArrayList<Individual>();
        for (int i = 0; i < caseNum; i++) {
            sampledIndiv.add(allIndiv.get(caseIndexes.getQuick((int) caseindexes[i])));
        }

        long[] controlindexes = new long[controlNum];
        RandomSampler.sample(controlNum, controlIndexes.size(), controlNum, 0, controlindexes, 0, new cern.jet.random.engine.MersenneTwister(new java.util.Date()));

        for (int i = 0; i < controlNum; i++) {
            sampledIndiv.add(allIndiv.get(controlIndexes.getQuick((int) controlindexes[i])));
        }

        return sampledIndiv;
    }

    public List<Individual> sampleWithoutReplacement(List<Individual> allIndiv, int caseNum) {
        //http://acs.lbl.gov/software/colt/api/cern/jet/random/sampling/RandomSampler.html#RandomSampler(long, long, long, cern.jet.random.engine.RandomEngine)       
        long[] caseindexes = new long[caseNum];
        RandomSampler.sample(caseNum, allIndiv.size(), caseNum, 0, caseindexes, 0, new cern.jet.random.engine.MersenneTwister(new java.util.Date()));
        List<Individual> sampledIndiv = new ArrayList<Individual>();
        for (int i = 0; i < caseNum; i++) {
            sampledIndiv.add(allIndiv.get((int) caseindexes[i]));
        }
        return sampledIndiv;
    }

    public double selectExtreemSample(List<Individual> allIndiv, double ratio) {
        int indSize = allIndiv.size();
        Collections.sort(allIndiv, new IndividualComparator());
        int selectSize = (int) (indSize * ratio);
        List<Individual> sampledIndiv = new ArrayList<Individual>();
        indSize--;
        DoubleArrayList orgList = new DoubleArrayList();
        DoubleArrayList selList = new DoubleArrayList();
        for (int i = 0; i < indSize; i++) {
            orgList.add(allIndiv.get(i).getMainTrait()[0]);
        }
        for (int i = 0; i < selectSize; i++) {
            sampledIndiv.add(allIndiv.get(i));
            sampledIndiv.add(allIndiv.get(indSize - i));
            selList.add(allIndiv.get(i).getMainTrait()[0]);
            selList.add(allIndiv.get(indSize - i).getMainTrait()[0]);
        }
        allIndiv.clear();
        allIndiv.addAll(sampledIndiv);
        double orgV = Descriptive.covariance(orgList, orgList);
        double selV = Descriptive.covariance(selList, selList);
        return selV / orgV;
    }

    public List<Individual> makeGenotypingErrors(List<Individual> allIndiv, int errorSite, double errorRate) {
        List<Individual> sampledIndiv = new ArrayList<Individual>();
        Iterator<Individual> iter = allIndiv.iterator();
        while (iter.hasNext()) {
            Individual indiv = iter.next();
            double r = randGenerator.nextUniform(0, 1);

            if (r <= errorRate) {
                //going to make errors
                Individual newIndiv = new Individual();
                newIndiv.setAffectedStatus(indiv.getAffectedStatus());
                newIndiv.markerGtySet = new StatusGtySet();
                newIndiv.markerGtySet.existence = indiv.markerGtySet.existence.copy();
                newIndiv.markerGtySet.maternalChrom = indiv.markerGtySet.maternalChrom.copy();
                newIndiv.markerGtySet.paternalChrom = indiv.markerGtySet.paternalChrom.copy();
                if (newIndiv.markerGtySet.existence.getQuick(errorSite)) {
                    r = randGenerator.nextUniform(0, 1);
                    if (r < 0.333) {
                        newIndiv.markerGtySet.maternalChrom.putQuick(errorSite, !newIndiv.markerGtySet.maternalChrom.getQuick(errorSite));
                        newIndiv.markerGtySet.paternalChrom.putQuick(errorSite, !newIndiv.markerGtySet.paternalChrom.getQuick(errorSite));
                    } else if (r < 0.6667) {
                        newIndiv.markerGtySet.maternalChrom.putQuick(errorSite, !newIndiv.markerGtySet.maternalChrom.getQuick(errorSite));
                    } else {
                        newIndiv.markerGtySet.paternalChrom.putQuick(errorSite, !newIndiv.markerGtySet.paternalChrom.getQuick(errorSite));
                    }
                }
                sampledIndiv.add(newIndiv);
            } else {
                sampledIndiv.add(indiv);
            }
        }

        return sampledIndiv;
    }

    public List<Individual> randomSplitSample(List<Individual> allIndiv, int caseNum, int controlNum) {
        Collections.shuffle(allIndiv);
        List<Individual> sampledIndiv = new ArrayList<Individual>();
        int sampleSize = allIndiv.size();
        for (int i = sampleSize - 1; i >= 0; i--) {
            Individual indiv = allIndiv.get(i);
            if (caseNum > 0 && indiv.getAffectedStatus() == 2) {
                sampledIndiv.add(allIndiv.remove(i));
                caseNum--;
            } else if (controlNum > 0 && indiv.getAffectedStatus() == 1) {
                sampledIndiv.add(allIndiv.remove(i));
                controlNum--;
            }
            if (caseNum == 0 && controlNum == 0) {
                break;
            }
        }

        if (caseNum != 0) {
            System.out.println("Warning! The population has less case than " + caseNum);
        }

        if (controlNum != 0) {
            System.out.println("Warning! The population has less control than " + controlNum);
        }
        return sampledIndiv;
    }

    public synchronized void randomSelectSample(List<Individual> allIndiv, int caseNum, int controlNum,
            List<Individual> selectedIndivs, List<Individual> remainderIndivs) {
        Collections.shuffle(allIndiv);
        int sampleSize = allIndiv.size();
        int i = sampleSize - 1;
        for (; i >= 0; i--) {
            Individual indiv = allIndiv.get(i);
            if (caseNum > 0 && indiv.getAffectedStatus() == 2) {
                selectedIndivs.add(indiv);
                caseNum--;
            } else if (controlNum > 0 && indiv.getAffectedStatus() == 1) {
                selectedIndivs.add(indiv);
                controlNum--;
            } else {
                remainderIndivs.add(indiv);
            }
            if (caseNum == 0 && controlNum == 0) {
                break;
            }
        }
        if (caseNum != 0) {
            System.out.println("Warning! The population has less case than " + caseNum);
        }
        if (controlNum != 0) {
            System.out.println("Warning! The population has less control than " + controlNum);
        }
        for (; i >= 0; i--) {
            remainderIndivs.add(allIndiv.get(i));
        }
    }

    public void permuteDiseaseStatus(List<Individual> allIndiv) {
        int N = allIndiv.size();
        int lessN = N - 1;
        for (int i = 0; i < lessN; i++) {
            int r = i + randGenerator.nextInt(0, lessN - i);   // between i and N-1
            int diseaes1 = allIndiv.get(i).getAffectedStatus();
            int diseaes2 = allIndiv.get(r).getAffectedStatus();
            allIndiv.get(i).setAffectedStatus(diseaes2);
            allIndiv.get(r).setAffectedStatus(diseaes1);
        }

    }

    public void randomAssignDiseaseStatus(List<Individual> allIndiv, double caseRatio) {
        int N = allIndiv.size();

        for (int i = 0; i < N; i++) {
            double r = randGenerator.nextUniform(0, 1);
            if (r < caseRatio) {
                allIndiv.get(i).setAffectedStatus(2);
            } else {
                allIndiv.get(i).setAffectedStatus(1);
            }
        }
    }

    public void randomAssignDiseaseStatusAdjustSex(List<Individual> allIndiv, double caseRatio) {
        int N = allIndiv.size();
        double PY2X1 = 0.1;
        double PY2 = 0.15, PX1 = 0.5;
        double PX1Y2 = PX1 / PY2 * PY2X1;
        double PX1Y1 = PX1 / (1 - PY2) * (1 - PY2X1);
        for (int i = 0; i < N; i++) {
            double r = randGenerator.nextUniform(0, 1);
            if (r < caseRatio) {
                allIndiv.get(i).setAffectedStatus(2);
                r = randGenerator.nextUniform(0, 1);
                if (r < PX1Y2) {
                    allIndiv.get(i).setGender(1);
                    allIndiv.get(i).getTraits().set(0, String.valueOf(1));
                } else {
                    allIndiv.get(i).setGender(2);
                    allIndiv.get(i).getTraits().set(0, String.valueOf(2));
                }
            } else {
                allIndiv.get(i).setAffectedStatus(1);
                r = randGenerator.nextUniform(0, 1);
                if (r < PX1Y1) {
                    allIndiv.get(i).setGender(1);
                    allIndiv.get(i).getTraits().set(0, String.valueOf(1));
                } else {
                    allIndiv.get(i).setGender(2);
                    allIndiv.get(i).getTraits().set(0, String.valueOf(2));
                }
            }
        }

    }

    public void randomAssignTwoBinaryTraits(List<Individual> allIndiv, double condiA, double condia, double caseRatio) {
        int N = allIndiv.size();

        for (int i = 0; i < N; i++) {
            double[] ts = new double[2];
            double r = randGenerator.nextUniform(0, 1);
            //A.
            if (r < caseRatio) {
                ts[0] = 2;
                r = randGenerator.nextUniform(0, 1);
                if (r < condiA) {
                    ts[1] = 2;
                } else {
                    ts[1] = 1;
                }

            } else {
                //a.
                ts[0] = 1;
                r = randGenerator.nextUniform(0, 1);
                if (r < condia) {
                    ts[1] = 2;
                } else {
                    ts[1] = 1;
                }
            }
            allIndiv.get(i).setMainTrait(ts);
        }

    }

    public void randomAssignQuantitativeTraits(List<Individual> allIndiv) {
        int N = allIndiv.size();
        for (int i = 0; i < N; i++) {
            double r = randGenerator.nextGaussian(0, 1);
            allIndiv.get(i).setMainTrait(new double[]{r});
        }
    }

    public void randomAssignQuantitativeTraitsBySex(List<Individual> allIndiv) {
        int N = allIndiv.size();
        double r = 0;
        for (int i = 0; i < N; i++) {
            if (randGenerator.nextUniform(0, 1) < 0.5) {
                r = randGenerator.nextGaussian(0, 1);
                allIndiv.get(i).setGender(1);
            } else {
                r = randGenerator.nextGaussian(0.2, 1);
                allIndiv.get(i).setGender(2);
            }
            allIndiv.get(i).setMainTrait(new double[]{r});
        }
    }

    public void randomAssignMultiQuantitativeTraits(List<Individual> allIndiv, double[] mean, DoubleMatrix2D covarianceMatrix) {
        final double PRECISION = 1.0e-8;
        int dim = mean.length;
        double[] tmpSample = new double[dim];
        //shame: I do not know how to use it
        NormalGen ng = new NormalGen(new MT19937(new WELL607()));
        MultinormalGen sg = new MultinormalCholeskyGen(ng, mean, covarianceMatrix);
        try {
            //CorrelatedRandomVectorGenerator sg = new CorrelatedRandomVectorGenerator(mean, covarianceMatrix, PRECISION, new GaussianRandomGenerator(new MersenneTwister()));
            int N = allIndiv.size();
            for (int i = 0; i < N; i++) {
                sg.nextPoint(tmpSample);
                double[] samples = new double[dim];
                System.arraycopy(tmpSample, 0, samples, 0, dim);
                allIndiv.get(i).setMainTrait(samples);
            }
        } catch (Exception ex) {
            ex.printStackTrace();
        }

    }

    public synchronized List<Individual> permuteNewDiseaseStatus(final List<Individual> allIndiv) {
        List<Individual> tmpInidivList = new ArrayList<Individual>();
        int N = allIndiv.size();
        int lessN = N - 1;

        //generate a new individual list which only includes disease information .
        for (int i = 0; i < N; i++) {
            Individual newIndiv = new Individual();
            newIndiv.setAffectedStatus(allIndiv.get(i).getAffectedStatus());
            tmpInidivList.add(newIndiv);
        }

        for (int i = 0; i < lessN; i++) {
            int r = i + randGenerator.nextInt(0, lessN - i);   // between i and N-1
            int diseaes1 = tmpInidivList.get(i).getAffectedStatus();
            int diseaes2 = tmpInidivList.get(r).getAffectedStatus();
            tmpInidivList.get(i).setAffectedStatus(diseaes2);
            tmpInidivList.get(r).setAffectedStatus(diseaes1);
        }

        return tmpInidivList;
    }
}
