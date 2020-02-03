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

import java.util.Comparator;

/**
 *
 * @author Miaoxin Li
 */
public class SNPPosiComparator implements Comparator {

    public int compare(Object arg0, Object arg1) {
        SNP obj1 = (SNP) arg0;
        SNP obj2 = (SNP) arg1;
        return (int) (obj1.getPhysicalPosition() - obj2.getPhysicalPosition());
    }
}
