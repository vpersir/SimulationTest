/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package org.cobi.genetsimulator.entity;

/**
 *
 * @author mxli
 */
public class AssociationSNP extends SNP{
    double logOddsRatio;
    boolean highRiskAllele=true;
    char oddsRatioAllle='X';
    public int orderInList=-1;
    double p;

    public double getP() {
        return p;
    }

    public void setP(double p) {
        this.p = p;
    }


    public AssociationSNP(String RSID, int chromosomeNum, int physicalPosition, double decodeGeneticMap) {
        super(RSID, chromosomeNum, physicalPosition, decodeGeneticMap);
    }

    public AssociationSNP(String RSID, int chromosomeNum, int physicalPosition) {
        super(RSID, chromosomeNum, physicalPosition);
    }

    public boolean isHighRiskAllele() {
        return highRiskAllele;
    }

    public void setHighRiskAllele(boolean highRiskAllele) {
        this.highRiskAllele = highRiskAllele;
    }


    public AssociationSNP() {
    }

    public double getLogOddsRatio() {
        return logOddsRatio;
    }

    public void setLogOddsRatio(double logOddsRatio) {
        this.logOddsRatio = logOddsRatio;
    }

    public char getOddsRatioAllle() {
        return oddsRatioAllle;
    }

    public void setOddsRatioAllle(char oddsRatioAllle) {
        this.oddsRatioAllle = oddsRatioAllle;
    }
     

}
