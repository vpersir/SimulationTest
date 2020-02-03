/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.cobi.genetsimulator.entity;

import java.io.Serializable;
import java.util.ArrayList;
import java.util.List;

/**
 *
 * @author Miaoxin Li
 */
public class Individual implements Cloneable, Serializable {

    private String familyID;
    private String individualID;
    private String momID;
    private String dadID;
    private int gender;
    private int affectedStatus;
    private int liability; //optional
    public StatusGtySet markerGtySet;
    public StatusGtySet traitGtySet;
    private List<String> traitValues;
    public double[] mainTrait;
    private String labelInChip;
    public StatusGtySet[] markerGtySets;
    private int chromNum = 24;

    public Individual() {
        traitValues = new ArrayList<String>();
        markerGtySets = new StatusGtySet[chromNum];
    }

    public Individual(String individualID) {
        this.individualID = individualID;
        traitValues = new ArrayList<String>();
        markerGtySets = new StatusGtySet[chromNum];
    }

    public double[] getMainTrait() {
        return mainTrait;
    }

    public void setMainTrait(double[] mainTrait) {
        this.mainTrait = mainTrait;
    }

    @Override
    public Object clone() throws CloneNotSupportedException {
        Individual o = null;
        o = (Individual) super.clone();
        o.traitValues = new ArrayList<String>();
        o.traitValues.addAll(this.traitValues);
        o.markerGtySet = new StatusGtySet();
        o.markerGtySet.existence = this.markerGtySet.existence.copy();
        o.markerGtySet.maternalChrom = this.markerGtySet.maternalChrom.copy();
        o.markerGtySet.paternalChrom = this.markerGtySet.paternalChrom.copy();

        return o;
    }

    public String getLabelInChip() {
        return labelInChip;
    }

    public StatusGtySet getMarkerGtySet() {
        return markerGtySet;
    }

    public void setMarkerGtySet(StatusGtySet markerGtySet) {
        this.markerGtySet = markerGtySet;
    }

    public StatusGtySet getTraitGtySet() {
        return traitGtySet;
    }

    public void setTraitGtySet(StatusGtySet traitGtySet) {
        this.traitGtySet = traitGtySet;
    }

    public void setLabelInChip(String labelInChip) {
        this.labelInChip = labelInChip;
    }

    /**
     * Get the value of traitValues
     *
     * @return the value of traitValues
     */
    public List<String> getTraits() {
        return traitValues;
    }

    /**
     * Set the value of traitValues
     *
     * @param traitValues new value of traitValues
     */
    public void setTraits(ArrayList<String> traits) {
        this.traitValues = traits;
    }

    /**
     * Set the value of a trait
     *
     * @param traitValues new value of trait
     */
    public void addTrait(String trait) {
        this.traitValues.add(trait);
    }

    public int getAffectedStatus() {
        return affectedStatus;
    }

    public void setAffectedStatus(int affectedStatus) {
        this.affectedStatus = affectedStatus;
    }

    public String getDadID() {
        return dadID;
    }

    public void setDadID(String dadID) {
        this.dadID = dadID;
    }

    public String getFamilyID() {
        return familyID;
    }

    public void setFamilyID(String familyID) {
        this.familyID = familyID;
    }

    public int getGender() {
        return gender;
    }

    public void setGender(int gender) {
        this.gender = gender;
    }

    public String getIndividualID() {
        return individualID;
    }

    public void setIndividualID(String individualID) {
        this.individualID = individualID;
    }

    public int getLiability() {
        return liability;
    }

    public void setLiability(int liability) {
        this.liability = liability;
    }

    public String getMomID() {
        return momID;
    }

    public void setMomID(String momID) {
        this.momID = momID;
    }
}
