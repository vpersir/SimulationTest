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

import java.io.BufferedInputStream;
import java.io.BufferedOutputStream;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.ObjectInputStream;
import java.io.ObjectOutputStream;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

/**
 *
 * @author Miaoxin Li
 */
public class ChromSNPGtySet {

    private String chromName = null;
    private List<AnnotSNP> snpInfoList = null;
    private Map<String, StatusGtySet> indiviGty = null;
    private int totalGtyBlockNum = 4;
    int curGtyBlockIndex = 0;

    public int getTotalGtyBlockNum() {
        return totalGtyBlockNum;
    }

    public void setTotalGtyBlockNum(int totalGtyBlockNum) {
        this.totalGtyBlockNum = totalGtyBlockNum;
    }

    public String getChromName() {
        return chromName;
    }

    public int getCurGtyBlockIndex() {
        return curGtyBlockIndex;
    }

    public void setCurGtyBlockIndex(int curGtyBlockIndex) {
        this.curGtyBlockIndex = curGtyBlockIndex;
    }

    public void setChromName(String chromName) {
        this.chromName = chromName;
    }

    public Map<String, StatusGtySet> getIndiviGty() {
        return indiviGty;
    }

    public void setIndiviGty(Map<String, StatusGtySet> linkPedGty) {
        this.indiviGty = linkPedGty;
    }

    public List<AnnotSNP> getSnpInfoList() {
        return snpInfoList;
    }

    public void setSnpInfoList(List<AnnotSNP> snpInfoMap) {
        this.snpInfoList = snpInfoMap;
    }

    public void saveAllObjectsAndClear(String path, int blockNum) throws Exception {
        File chekPath = new File(path);
        if (!chekPath.exists()) {
            chekPath.mkdir();
        }
        FileOutputStream objFOut = new FileOutputStream(path + File.separator + chromName + ".snp");
        BufferedOutputStream objOBfs = new BufferedOutputStream(objFOut);
        ObjectOutputStream localObjOut = new ObjectOutputStream(objOBfs);
        totalGtyBlockNum = blockNum;
        //Serialization starts here.
        localObjOut.writeObject(this.chromName);
        localObjOut.writeObject(this.snpInfoList);
        this.snpInfoList.clear();
        int cutOff = indiviGty.size() / totalGtyBlockNum;
        if (cutOff == 0) {
            localObjOut.writeInt(1);
        } else {
            localObjOut.writeInt(totalGtyBlockNum);
        }
        localObjOut.flush();
        localObjOut.close();
        objOBfs.flush();
        objOBfs.close();
        objFOut.close();

        int index = 0;
        int curGtyBlockNum = 0;
        Map<String, StatusGtySet> tmpMap = new HashMap<String, StatusGtySet>();
        totalGtyBlockNum--;
        for (Map.Entry<String, StatusGtySet> m : indiviGty.entrySet()) {
            tmpMap.put(m.getKey(), m.getValue());
            index++;
            if (index == cutOff && curGtyBlockNum < totalGtyBlockNum) {
                objFOut = new FileOutputStream(path + File.separator + chromName + ".gty." + (curGtyBlockNum));
                objOBfs = new BufferedOutputStream(objFOut);
                localObjOut = new ObjectOutputStream(objOBfs);
                localObjOut.writeObject(tmpMap);
                localObjOut.flush();
                localObjOut.close();
                objOBfs.flush();
                objOBfs.close();
                objFOut.close();
                tmpMap = new HashMap<String, StatusGtySet>();
                index = 0;
                curGtyBlockNum++;
            }
        }
        totalGtyBlockNum++;
        objFOut = new FileOutputStream(path + File.separator + chromName + ".gty." + (curGtyBlockNum));
        objOBfs = new BufferedOutputStream(objFOut);
        localObjOut = new ObjectOutputStream(objOBfs);
        localObjOut.writeObject(tmpMap);
        localObjOut.flush();
        localObjOut.close();
        objOBfs.flush();
        objOBfs.close();
        objFOut.close();
        tmpMap.clear();
        this.indiviGty.clear();

    }

    public void saveAllBufObjectsAndClear(String path, int bufNum) throws Exception {
        FileOutputStream objFOut = new FileOutputStream(path + File.separator + chromName + ".snp.buf" + bufNum);
        BufferedOutputStream objOBfs = new BufferedOutputStream(objFOut);
        ObjectOutputStream localObjOut = new ObjectOutputStream(objOBfs);

        //Serialization starts here.
        localObjOut.writeObject(this.chromName);
        localObjOut.writeObject(this.snpInfoList);
        this.snpInfoList.clear();
        //1 block for all genotypes
        localObjOut.writeInt(1);
        localObjOut.flush();
        localObjOut.close();
        objOBfs.flush();
        objOBfs.close();
        objFOut.close();

        objFOut = new FileOutputStream(path + File.separator + chromName + ".gty.buf" + bufNum);
        objOBfs = new BufferedOutputStream(objFOut);
        localObjOut = new ObjectOutputStream(objOBfs);
        localObjOut.writeObject(this.indiviGty);
        localObjOut.flush();
        localObjOut.close();
        objOBfs.flush();
        objOBfs.close();
        objFOut.close();
        this.indiviGty.clear();
    }

    public void readAllBufObjectsAndDelete(String path, int bufNum) throws Exception {
        File fileName = new File(path + File.separator + chromName + ".snp.buf" + bufNum);
        FileInputStream objFIn = new FileInputStream(fileName);
        BufferedInputStream objIBfs = new BufferedInputStream(objFIn);
        ObjectInputStream localObjIn = new ObjectInputStream(objIBfs);
        this.chromName = (String) localObjIn.readObject();
        this.snpInfoList = (List<AnnotSNP>) localObjIn.readObject();
        this.totalGtyBlockNum = localObjIn.readInt();
        curGtyBlockIndex = 0;
        localObjIn.close();
        objIBfs.close();
        objFIn.close();
        fileName.delete();

        fileName = new File(path + File.separator + chromName + ".gty.buf" + bufNum);
        objFIn = new FileInputStream(fileName);
        objIBfs = new BufferedInputStream(objFIn);
        localObjIn = new ObjectInputStream(objIBfs);
        this.indiviGty = (Map<String, StatusGtySet>) localObjIn.readObject();

        localObjIn.close();
        objIBfs.close();
        objFIn.close();
        fileName.delete();
    }

    public void readAllObjects(String path) throws Exception {
        File checkFile = new File(path + File.separator + chromName + ".snp");
        if (!checkFile.exists()) {
            this.snpInfoList = new ArrayList<AnnotSNP>();
            this.indiviGty = new HashMap<String, StatusGtySet>();
            return;
        }
        FileInputStream objFIn = new FileInputStream(checkFile);
        BufferedInputStream objIBfs = new BufferedInputStream(objFIn);
        ObjectInputStream localObjIn = new ObjectInputStream(objIBfs);
        this.chromName = (String) localObjIn.readObject();
        this.snpInfoList = (List<AnnotSNP>) localObjIn.readObject();
        this.totalGtyBlockNum = localObjIn.readInt();
        curGtyBlockIndex = 0;
        localObjIn.close();
        objIBfs.close();
        objFIn.close();

        this.indiviGty = new HashMap<String, StatusGtySet>();

        for (int i = 0; i < totalGtyBlockNum; i++) {
            objFIn = new FileInputStream(path + File.separator + chromName + ".gty." + curGtyBlockIndex);
            objIBfs = new BufferedInputStream(objFIn);
            localObjIn = new ObjectInputStream(objIBfs);

            Map<String, StatusGtySet> tmpMap = (Map<String, StatusGtySet>) localObjIn.readObject();
            if (tmpMap.size() == 0) {
                break;
            }
            curGtyBlockIndex++;
            localObjIn.close();
            objIBfs.close();
            objFIn.close();
            this.indiviGty.putAll(tmpMap);
        }
    }

    public void saveChromSNPInforAndClear(String path) throws Exception {
        FileOutputStream objFOut = new FileOutputStream(path + File.separator + chromName + ".snp");
        BufferedOutputStream objOBfs = new BufferedOutputStream(objFOut);
        ObjectOutputStream localObjOut = new ObjectOutputStream(objOBfs);
        //Serialization starts here.
        localObjOut.writeObject(this.chromName);
        localObjOut.writeObject(this.snpInfoList);
        this.snpInfoList.clear();
        localObjOut.writeInt(totalGtyBlockNum);
        localObjOut.flush();
        localObjOut.close();
        objOBfs.flush();
        objOBfs.close();
        objFOut.close();
    }

    public void saveGtyBlock(String path) throws Exception {
        FileOutputStream objFOut = new FileOutputStream(path + File.separator + chromName + ".gty." + (curGtyBlockIndex - 1));
        BufferedOutputStream objOBfs = new BufferedOutputStream(objFOut);
        ObjectOutputStream localObjOut = new ObjectOutputStream(objOBfs);
        localObjOut.writeObject(this.indiviGty);
        localObjOut.flush();
        localObjOut.close();
        objOBfs.flush();
        objOBfs.close();
        objFOut.close();
    }

    //readChromSNPInfo, hasMoreGtyBlock and readNextGtyBlock work together to read data part by part
    public void readChromSNPInfo(String path) throws Exception {
        FileInputStream objFIn = new FileInputStream(path + File.separator + chromName + ".snp");
        BufferedInputStream objIBfs = new BufferedInputStream(objFIn);
        ObjectInputStream localObjIn = new ObjectInputStream(objIBfs);
        this.chromName = (String) localObjIn.readObject();
        this.snpInfoList = (List<AnnotSNP>) localObjIn.readObject();
        this.totalGtyBlockNum = localObjIn.readInt();
        curGtyBlockIndex = 0;
        localObjIn.close();
        objIBfs.close();
        objFIn.close();
    }

    public boolean hasMoreGtyBlock() {
        return (curGtyBlockIndex < totalGtyBlockNum);
    }

    public void readNextGtyBlock(String path) throws Exception {
        FileInputStream objFIn = new FileInputStream(path + File.separator + chromName + ".gty." + curGtyBlockIndex);
        BufferedInputStream objIBfs = new BufferedInputStream(objFIn);
        ObjectInputStream localObjIn = new ObjectInputStream(objIBfs);
        this.indiviGty = (Map<String, StatusGtySet>) localObjIn.readObject();
        curGtyBlockIndex++;
        localObjIn.close();
        objIBfs.close();
        objFIn.close();
    }
}
