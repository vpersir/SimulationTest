/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package org.cobi.genetsimulator.entity;

import java.io.Serializable;
import org.cobi.eqtlsimulator.Constants;
import static org.cobi.genetsimulator.entity.PlinkDataset.MISSING_ALLELE_NAME;
import static org.cobi.genetsimulator.entity.PlinkDataset.MISSING_STRAND_NAME;
 

/**
 *
 * @author mxli
 */
public class AnnotSNP extends SNP implements Serializable,Constants {
    //in mapfile it indicates the first genotype positions in pedigree file, 0, 1,2,3 ...
    public int probeSetID;
    private String flankingSequence;
    private char aAllele=MISSING_ALLELE_NAME;
    private char bAllele=MISSING_ALLELE_NAME;
    private double aAlleleFreq;
  
    private char strand=MISSING_STRAND_NAME;
    //indicate the first genotype positions in the bitvector, 0, 2, 4, 6, 8 ... 
    public int order;

    public int index;

    public char getaAllele() {
        return aAllele;
    }

    public void setaAllele(char aAllele) {
        this.aAllele = aAllele;
    }

    public double getaAlleleFreq() {
        return aAlleleFreq;
    }

    public void setaAlleleFreq(double aAlleleFreq) {
        this.aAlleleFreq = aAlleleFreq;
    }


    public AnnotSNP(String RSID, int chromosomeNum, int physicalPosition, double decodeGeneticMap, int probeSetID, String flankingSequence, char aAllele, char bAllele, double aAlleleFreq, double bAlleleFreq, int order, char strand) {
        super(RSID, chromosomeNum, physicalPosition, decodeGeneticMap);
        this.probeSetID = probeSetID;
        this.flankingSequence = flankingSequence;
        this.aAllele = aAllele;
        this.bAllele = bAllele;
        this.aAlleleFreq = aAlleleFreq;
        this.order = order;
        this.strand = strand;
    }

    public char getStrand() {
        return strand;
    }

    public void setStrand(char strand) {
        this.strand = strand;
    }

    public AnnotSNP(String RSID, int chromosomeNum, int physicalPosition, double decodeGeneticMap) {
        super(RSID, chromosomeNum, physicalPosition, decodeGeneticMap);
    }

    public AnnotSNP(String RSID, int chromosomeNum, int physicalPosition, double decodeGeneticMap, char aAllele, char bAllele) {
        super(RSID, chromosomeNum, physicalPosition, decodeGeneticMap);
        this.aAllele = aAllele;
        this.bAllele = bAllele;
    }

    public AnnotSNP() {
    }

    public AnnotSNP(String RSID, int chromosomeNum, int physicalPosition, double decodeGeneticMap, int probeSetID, String flankingSequence, char aAllele, char bAllele, double aAlleleFreq, double bAlleleFreq) {
        super(RSID, chromosomeNum, physicalPosition, decodeGeneticMap);
        this.probeSetID = probeSetID;
        this.flankingSequence = flankingSequence;
        this.aAllele = aAllele;
        this.bAllele = bAllele;
        this.aAlleleFreq = aAlleleFreq;
    }

    /**
     * Get the value of aAlleleFreq
     *
     * @return the value of aAlleleFreq
     */
    public double getAAlleleFreq() {
        return aAlleleFreq;
    }

    /**
     * Set the value of aAlleleFreq
     *
     * @param aAlleleFreq new value of aAlleleFreq
     */
    public void setAAlleleFreq(double aAlleleFreq) {
        this.aAlleleFreq = aAlleleFreq;
    }

    /**
     * Get the value of bAllele
     *
     * @return the value of bAllele
     */
    public char getBAllele() {
        return bAllele;
    }

    /**
     * Set the value of bAllele
     *
     * @param bAllele new value of bAllele
     */
    public void setBAllele(char bAllele) {
        this.bAllele = bAllele;
    }

    /**
     * Get the value of Aallele
     *
     * @return the value of Aallele
     */
    public char getAAllele() {
        return aAllele;
    }

    /**
     * Set the value of Aallele
     *
     * @param Aallele new value of Aallele
     */
    public void setAAllele(char Aallele) {
        this.aAllele = Aallele;
    }

    /**
     * Get the value of flankingSequence
     *
     * @return the value of flankingSequence
     */
    public String getFlankingSequence() {
        return flankingSequence;
    }

    /**
     * Set the value of flankingSequence
     *
     * @param flankingSequence new value of flankingSequence
     */
    public void setFlankingSequence(String flankingSequence) {
        this.flankingSequence = flankingSequence;
    }

    /**
     * Get the value of probeSetID
     *
     * @return the value of probeSetID
     */
    public int getProbeSetID() {
        return probeSetID;
    }

    /**
     * Set the value of probeSetID
     *
     * @param probeSetID new value of probeSetID
     */
    public void setProbeSetID(int probeSetID) {
        this.probeSetID = probeSetID;
    }
}