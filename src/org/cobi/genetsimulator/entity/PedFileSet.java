/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.cobi.genetsimulator.entity;

/**
 *
 * @author mxli
 */
public class PedFileSet {

    protected String pedigreeFileName;
    protected String mapFileName;

    public PedFileSet(String pedigreeFileName, String mapFileName) {
        this.pedigreeFileName = pedigreeFileName;
        this.mapFileName = mapFileName;
    }

    public PedFileSet() {
    }

    public String getMapFileName() {
        return mapFileName;
    }

    public void setMapFileName(String mapFileName) {
        this.mapFileName = mapFileName;
    }

    public String getPedigreeFileName() {
        return pedigreeFileName;
    }

    public void setPedigreeFileName(String pedigreeFileName) {
        this.pedigreeFileName = pedigreeFileName;
    }
 
    
    
}
