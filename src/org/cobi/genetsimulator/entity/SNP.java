// (c) 2008-2009 Miaoxin Li
// This file is distributed as part of the IGG source code package
// and may not be redistributed in any form, without prior written
// permission from the author. Permission is granted for you to
// modify this file for your own personal use, but modified versions
// must retain this copyright notice and must not be distributed.

// Permission is granted for you to use this file to compile IGG.

// All computer programs have bugs. Use this file at your own risk.
// Saturday, January 17, 2009
/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.cobi.genetsimulator.entity;

import java.io.Serializable;

/**
 *
 * @author mxli
 */
public class SNP implements Serializable{

    protected String RSID;
    //X,Y, and Mitochondria are 23, 24,25
    protected int chromosomeNum;    //unknown postition is -9
    public int physicalPosition;
    //unknown map is -9.0
    protected double decodeGeneticMap;

    public SNP() {
    }

    public SNP(String RSID, int chromosomeNum, int physicalPosition) {
        this.RSID = RSID;
        this.chromosomeNum = chromosomeNum;
        this.physicalPosition = physicalPosition;
    }

    
    public SNP(String RSID, int chromosomeNum, int physicalPosition, double decodeGeneticMap) {
        this.RSID = RSID;
        this.chromosomeNum = chromosomeNum;
        this.physicalPosition = physicalPosition;
        this.decodeGeneticMap = decodeGeneticMap;
    }

    /**
     * Get the value of decodeGeneticMap
     *
     * @return the value of decodeGeneticMap
     */
    public double getDecodeGeneticMap() {
        return decodeGeneticMap;
    }

    /**
     * Set the value of decodeGeneticMap
     *
     * @param decodeGeneticMap new value of decodeGeneticMap
     */
    public void setDecodeGeneticMap(double decodeGeneticMap) {
        this.decodeGeneticMap = decodeGeneticMap;
    }

    /**
     * Get the value of physicalPosition
     *
     * @return the value of physicalPosition
     */
    public int getPhysicalPosition() {
        return physicalPosition;
    }

    /**
     * Set the value of physicalPosition
     *
     * @param physicalPosition new value of physicalPosition
     */
    public void setPhysicalPosition(int physicalPosition) {
        this.physicalPosition = physicalPosition;
    }

    /**
     * Get the value of RSID
     *
    private int chromosomeNum;
    
    /**
     * Get the value of chromosomeNum
     *
     * @return the value of chromosomeNum
     */
    public int getChromosomeNum() {
        return chromosomeNum;
    }

    /**
     * Set the value of chromosomeNum
     *
     * @param chromosomeNum new value of chromosomeNum
     */
    public void setChromosomeNum(int chromosomeNum) {
        this.chromosomeNum = chromosomeNum;
    }
    /*
    @return
    the value of RSID   
     */

    public String getRSID() {
        return RSID;
    }

    /**
     * Set the value of RSID
     *
     * @param RSID new value of RSID
     */
    public void setRSID(String RSID) {
        this.RSID = RSID;
    }
}
