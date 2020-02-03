/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.cobi.util.thread;

import java.util.ArrayList;
import java.util.List;

/**
 *
 * @author mxli
 */
public class RWResource {

    /** 
     * When readerNumber == 0, there is no one reading or writing. 
     * When readerNumber > 0, readerNumber means number of readers. 
     * When readerNumber < 0, it means that some writer is writing. */
    private int readerNumber = 0;
    /** 
     * the shared resource for writing or reading 
     */
    private List buffer = null;

    public RWResource() {
        buffer = new ArrayList(512);
        readerNumber = 0;
    }

    /** 
     * get buffer for reading. 
     * should be called before reading 
     * @return the buffer 
     */
    public synchronized List getBufferForReading() {
// if some writer is writing, wait until no writer is writing 
        while (readerNumber < 0) {
            try {
                this.wait();
            } catch (InterruptedException e) {
                e.printStackTrace();
            }
        }
// when readerNumber >= 0 
        readerNumber++;
        return buffer;
    }

    /** 
     * should be called after reading 
     */
    public synchronized void finishReading() {
        readerNumber--;
        if (readerNumber == 0) {
            this.notifyAll(); // notify possible waiting writers 
        }
    }

    /** 
     * get buffer for writing. 
     * should be called before writing. 
     * @return the buffer 
     */
    public synchronized List getBufferForWriting() {
// if some writer is writing or some reader is reading, wait until no one is writing or reading 
        while (readerNumber != 0) {
            try {
                this.wait();
            } catch (InterruptedException e) {
                e.printStackTrace();
            }
        }
// when readerNumber == 0 
        readerNumber--; // now readderNumber == -1. 
        return buffer;
    }

    /** 
     * should be called after writing */
    public synchronized void finishWriting() {
        readerNumber++; // readerNumber = -1 + 1 = 0; 
// readerNumber must be 0 at this point 
        this.notifyAll(); // notify possible waiting writers or waiting readers 
    }

    public static void main(String[] args) {
// init 
        RWResource resource = new RWResource();
// 
        RunnableReader reader = new RunnableReader();
        reader.setRWResource(resource);
// 
        RunnableWriter writer = new RunnableWriter();
        writer.setRWResource(resource);
        int writerNumber = 5; // 

        int readerNumber = 10; // 1 
// start writers  

        for (int i = 0; i < writerNumber; i++) {
            Thread thread = new Thread(writer, "writer" + (i + 1));
            thread.start();
        }
// give writers enough time to think and write articles 
        try {
            Thread.sleep(1000);
        } catch (InterruptedException e) {
            e.printStackTrace();
        }
 
        for (int i = 0; i < readerNumber; i++) {
            Thread thread = new Thread(reader, "reader" + (i + 1));
            thread.start();
        }
    }
}
 