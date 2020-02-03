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

import cern.colt.bitvector.BitVector;
import java.io.Serializable;

/**
 *
 * @author Miaoxin Li
 */
public class StatusGtySet implements Cloneable, Serializable {

    public BitVector paternalChrom;
    public BitVector maternalChrom;
    //0 indicating missing
    public BitVector existence;

    @Override
    protected Object clone() throws CloneNotSupportedException {
        StatusGtySet o = (StatusGtySet) super.clone();
        return o;
    }

    public StatusGtySet() {
    }

    public StatusGtySet(int size) {
        existence = new BitVector(size);
        paternalChrom = new BitVector(size);
        maternalChrom = new BitVector(size);
    }
}
