/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.cobi.genetsimulator.controller;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileReader;
import java.io.FileWriter;
import java.util.HashMap;
import java.util.Map;

/**
 *
 * @author mxli
 */
public class GeneInforProcessor {

    public Map<String, int[]> readLinkedPPIPairs(String genePosFile, String PPIFile, int extendLen) throws Exception {
        //read gene position information
        Map<String, String[]> geneInfoMap = new HashMap<String, String[]>();
        BufferedReader br = new BufferedReader(new FileReader(genePosFile));
        String line = null;
        //skip the head
        br.readLine();
        while ((line = br.readLine()) != null) {
            //line = line.trim();
            //System.out.println(line);
            if (line.trim().length() == 0) {
                continue;
            }
            String[] cells = line.split("\t", -1);
            cells[3] = String.valueOf(Integer.parseInt(cells[3]) - extendLen);
            cells[4] = String.valueOf(Integer.parseInt(cells[4]) + extendLen);
            geneInfoMap.put(cells[1], cells);
        }
        br.close();
        int difChromNum = 0;
        int sameChromNum = 0;
        int closeNum = 0;
        int closeDis = 100000;
        Map<String, int[]> genePositions = new HashMap<String, int[]>();
        br = new BufferedReader(new FileReader(PPIFile));
        while ((line = br.readLine()) != null) {
            //line = line.trim();
            //System.out.println(line);
            if (line.trim().length() == 0) {
                continue;
            }
            String[] cells = line.split("\t", -1);
            String geneA = cells[4];
            String geneB = cells[5];
            String[] geneACells = geneInfoMap.get(geneA);
            String[] geneBCells = geneInfoMap.get(geneB);
            if (geneACells == null) {
                // System.out.println(geneA + " has not position information!");
            }
            if (geneBCells == null) {
                // System.out.println(geneB + " has not position information!");
            }

            if (geneACells != null && geneBCells != null) {
                if (geneACells[2].equals(geneBCells[2])) {
                    sameChromNum++;
                    if (Integer.parseInt(geneACells[3]) < Integer.parseInt(geneBCells[3])) {
                        if (Math.abs(Integer.parseInt(geneBCells[3]) - Integer.parseInt(geneACells[4])) <= closeDis) {
                            closeNum++;
                            genePositions.put(geneA + ',' + geneB + ',' + geneACells[2],
                                    new int[]{Integer.parseInt(geneACells[3]), Integer.parseInt(geneACells[4]),
                                        Integer.parseInt(geneBCells[3]), Integer.parseInt(geneBCells[4])});
                        }
                    } else if (Math.abs(Integer.parseInt(geneACells[3]) - Integer.parseInt(geneBCells[4])) <= closeDis) {
                        closeNum++;
                        genePositions.put(geneB + ',' + geneA + ',' + geneACells[2],
                                new int[]{Integer.parseInt(geneBCells[3]), Integer.parseInt(geneBCells[4]),
                                    Integer.parseInt(geneACells[3]), Integer.parseInt(geneACells[4])});
                    }
                } else {
                    difChromNum++;
                }
            }
        }
        br.close();
        System.out.println("difChromNum " + difChromNum);
        System.out.println("sameChromNum " + sameChromNum);
        System.out.println("closeNum " + closeNum);
        return genePositions;
    }

    public Map<String, int[]> readConsecuiveGenePairs(String genePosFile, int extendLen, String chrName) throws Exception {
        //read gene position information

        BufferedReader br = new BufferedReader(new FileReader(genePosFile));
        String line = null;
        //skip the head
        br.readLine();
        line = br.readLine();
        Map<String, int[]> genePositions = new HashMap<String, int[]>();
        String[] cells0 = line.split("\t", -1);
        cells0[3] = String.valueOf(Integer.parseInt(cells0[3]) - extendLen);
        cells0[4] = String.valueOf(Integer.parseInt(cells0[4]) + extendLen);
        int closeDis = 1000000;
        line = br.readLine();
        do {
            //line = line.trim();
            //System.out.println(line);
            if (line.trim().length() == 0) {
                continue;
            }
            String[] cells = line.split("\t", -1);
            cells[3] = String.valueOf(Integer.parseInt(cells[3]) - extendLen);
            cells[4] = String.valueOf(Integer.parseInt(cells[4]) + extendLen);
            if (!cells[2].equals(chrName)) {
                break;
            }
            if (Integer.parseInt(cells0[3]) < Integer.parseInt(cells[3])) {
                if (Math.abs(Integer.parseInt(cells[3]) - Integer.parseInt(cells0[4])) <= closeDis) {

                    genePositions.put(cells0[1] + ',' + cells[1] + ',' + cells0[2],
                            new int[]{Integer.parseInt(cells0[3]), Integer.parseInt(cells0[4]),
                                Integer.parseInt(cells[3]), Integer.parseInt(cells[4])});
                }
            } else if (Math.abs(Integer.parseInt(cells0[3]) - Integer.parseInt(cells[4])) <= closeDis) {
                genePositions.put(cells[1] + ',' + cells0[1] + ',' + cells0[2],
                        new int[]{Integer.parseInt(cells[3]), Integer.parseInt(cells[4]),
                            Integer.parseInt(cells0[3]), Integer.parseInt(cells0[4])});
            }
            System.arraycopy(cells, 0, cells0, 0, cells.length);
        } while ((line = br.readLine()) != null);
        br.close();
        return genePositions;
    }

    public Map<String, int[]> readGenePos(String genePosFile, int extendLen, String chrName) throws Exception {
        //read gene position information
        BufferedReader br = new BufferedReader(new FileReader(genePosFile));
        String line = null;
        //skip the head
        br.readLine();
        line = br.readLine();
        Map<String, int[]> genePositions = new HashMap<String, int[]>();

        int closeDis = 1000000;
        line = br.readLine();
        do {
            //line = line.trim();
            //System.out.println(line);
            if (line.trim().length() == 0) {
                continue;
            }
            String[] cells = line.split("\t", -1);
            cells[3] = String.valueOf(Integer.parseInt(cells[3]) - extendLen);
            cells[4] = String.valueOf(Integer.parseInt(cells[4]) + extendLen);
            if (!cells[2].equals(chrName)) {
                break;
            }
            genePositions.put(cells[1] + "," + chrName, new int[]{Integer.parseInt(cells[3]), Integer.parseInt(cells[4])});
        } while ((line = br.readLine()) != null);
        br.close();
        return genePositions;
    }

    public void extractGenePositions(String readPath, String writePath) throws
            Exception {
        BufferedReader br = null;
        try {
            System.out.println("Parsing " + readPath + ". Please wait ...");

            br = new BufferedReader(new FileReader(readPath));
            BufferedWriter bwSeq_Gene_Map = new BufferedWriter(new FileWriter(writePath));
            String line = null;
            StringBuffer title = new StringBuffer();
            StringBuffer description = new StringBuffer();
            String[] linecells = null;
            if ((line = br.readLine()) != null) //skip the first line
            {
                while ((line = br.readLine()) != null) {
                    linecells = line.split("\t");
                    if ((linecells[11].equals("GENE")) && (linecells[12].startsWith("GRCh37"))) {
                        bwSeq_Gene_Map.write(linecells[10].substring(7)); //entreze gene id
                        bwSeq_Gene_Map.write('\t');
                        bwSeq_Gene_Map.write(linecells[9]); //entreze symbol
                        bwSeq_Gene_Map.write('\t');
                        bwSeq_Gene_Map.write(linecells[1].split("[|]")[0]);  //chromosome
                        bwSeq_Gene_Map.write('\t');
                        bwSeq_Gene_Map.write(linecells[2]); //chr_start
                        bwSeq_Gene_Map.write('\t');
                        bwSeq_Gene_Map.write(linecells[3]); //chr_stop
                        bwSeq_Gene_Map.write('\t');
                        bwSeq_Gene_Map.write(linecells[4]); //chr_orient
                        bwSeq_Gene_Map.write("\n");
                    }
                }
            }
            bwSeq_Gene_Map.close();
        } catch (Exception e) {
            e.printStackTrace();

        } finally {
            try {
                if (br != null) {
                    br.close();
                }
            } catch (Exception e) {
                e.printStackTrace();
            }
        }
    }
}
