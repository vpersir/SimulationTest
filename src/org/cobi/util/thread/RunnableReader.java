/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.cobi.util.thread;

import java.util.Iterator;
import java.util.List;

/**
 *
 * @author mxli
 */
public class RunnableReader implements Runnable {

    private RWResource resource = null;

    public RunnableReader() {
    }

    /** 
     * must be called before start running 
     * @param theResource 
     */
    public void setRWResource(RWResource theResource) {
        resource = theResource;
    }

    public void run() {
        while (true) {
// get the reader's name 
            String readerName = "[" + Thread.currentThread().getName() + "] ";
// first, get buffer for reading 
            List buffer = resource.getBufferForReading();
// reading 
            for (Iterator iterator = buffer.iterator(); iterator.hasNext();) {
                System.out.println(readerName + iterator.next());
            }
            int articleNumber = buffer.size();
            int thinkingTime = articleNumber * 1000;
            for (int i = 0; i < thinkingTime; i++) {
// thingking hard when reading 
                }
// finish reading 
            resource.finishReading();
// rest 
            try {
                Thread.sleep(articleNumber * 50);
            } catch (InterruptedException e) {
                e.printStackTrace();
            }
        }
    }
}