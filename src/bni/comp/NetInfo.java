/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package bni.comp;

import bni.Utils;
import mod.jmut.core.comp.NetData;

/**
 *
 * @author colin_PC2
 */
public class NetInfo {
    public NetData netD;
    
    //Scoring
    public double avg_score;
    public double[] g_scores;
        // the average score for each gene
    public double fitness;
    
    //Constraint
    public int[][] posPaths;
    public int[][] negPaths;
    
    public NetInfo(NetData netD, int[][] _posPaths, int[][] _negPaths) {
        int numGenes = netD.nodes.size();
        this.netD = netD;
        
        this.g_scores = new double[numGenes];
        java.util.Arrays.fill(this.g_scores, Double.MAX_VALUE);
        this.avg_score = Double.MAX_VALUE;
        this.fitness = -1;
                
        posPaths = new int[numGenes][numGenes];
        negPaths = new int[numGenes][numGenes];
        
        Utils.copy(_posPaths, posPaths);
        Utils.copy(_negPaths, negPaths);
    }
    
    public void average(int noExperiments) {
        double score = Utils.sum(this.g_scores);
        this.avg_score = score / (this.g_scores.length * noExperiments);
        
        this.fitness = 1 - this.avg_score;
    }
    
    public void add(double[] scores) {
        for(int i = 0; i < this.g_scores.length; i ++) {
            this.g_scores[i] = this.g_scores[i] + scores[i];
        }
    }
    
    public void free() {
        this.netD = null;
        this.g_scores = null;
        this.posPaths = null;
        this.negPaths = null;
    }
}
