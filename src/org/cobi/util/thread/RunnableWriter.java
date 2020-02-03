/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.cobi.util.thread;

import java.util.List;

/**
 *
 * @author mxli
 */
public class RunnableWriter implements Runnable {

    private RWResource resource = null;
    private int articleNumber = 0;

    public RunnableWriter() {
        articleNumber = 0;
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
// get the writer's name 
            String writerName = "[" + Thread.currentThread().getName() + "] ";
// first, get buffer for reading 
            List buffer = resource.getBufferForWriting();
            int nWritten = 3; // write 4 articles one time 

            for (int n = 0; n < nWritten; n++) {
// writing 
                articleNumber++;
                String articleName = "article" + articleNumber;
                buffer.add(articleName);
                System.out.println(writerName + articleName);
                int thinkingTime = 10000;
                for (int i = 0; i < thinkingTime; i++) {// thingking hard when writing 
                    }
            } // finish writing 

            resource.finishWriting();
// rest 
            try {
                Thread.sleep(500);
            } catch (InterruptedException e) {
                e.printStackTrace();
            }
        }
    }
}
