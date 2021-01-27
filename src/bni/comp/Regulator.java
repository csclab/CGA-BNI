/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package bni.comp;

/**
 *
 * @author colin_PC2
 */
public class Regulator implements Comparable<Regulator>{
    public static final int MODE_ADD = 0;
    public static final int MODE_DEL = 1;
    public static final int MODE_CHG = 2;
    
    public static final int RANK_DIRECT_POSITIVE = 1;
    public static final int RANK_DIRECT_NEGATIVE = 2;
    public static final int RANK_INDIRECT_POSITIVE = 3; // + > -
    public static final int RANK_INDIRECT_NEGATIVE = 4; // + < -
    public static final int RANK_INDIRECT_NEUTRAL = 5;  // + = -
    public static final int RANK_NOLINK = 6;
    public static final int RANK_UNKNOWN = 7;
    
    public Integer regulator;
    public Integer interaction;
    public int rank;
        /*
        *   1: differential expression value >= DIFF_THRESHOLD
        *   2: pearson correlation coefficient >= CORR_THRESHOLD
        * 
        * 
        * 
        */
    
    public Regulator(int srcGene, int type, int rank) {
        this.regulator = srcGene;
        this.interaction = type;
        this.rank = rank;
    }
    
    public String getLink() {        
        String s = " -> ";
        if(interaction == -1) s = " -| ";
        
        return String.valueOf(regulator) + s;
    }
    
    @Override
    public int compareTo(Regulator anotherLink) {
        int v = this.regulator.compareTo(anotherLink.regulator);
//        if (v != 0) return v;
//
//        v = this.interaction.compareTo(anotherLink.interaction);
        return v;
    }
        
}