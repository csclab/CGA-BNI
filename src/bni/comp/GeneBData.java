/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package bni.comp;

import bni.InferBN;
import bni.Utils;
import java.io.FileNotFoundException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.TreeSet;
import mod.jmut.core.comp.LogicTable;
import mod.jmut.core.comp.Node;

/**
 *
 * @author colin_PC2
 */
public class GeneBData {        
    public int numGenes;        
    
    boolean[][] data = null;
    int[] expGenes = null;
        //List of genes which are mutated
    int[] wildTypes = null;
    
    private HashMap<Integer, TreeSet<Regulator>> regulators =null;
        //key: gene id, from [1...numGenes]
        //value: regulators genes of the 'key' gene
    public double[][] corrs;
        //correlations between columns of genes
    
    public int[][] C;
        //Constraint matrix of rankings: [regulator][target]
    
    public boolean[] fixedGenes;
    public boolean finalSourceDest;//only perturb interaction sign
    
    public GeneBData(int numGenes) {
        this.numGenes = numGenes;                
        this.finalSourceDest = false;
                
        this.regulators = new HashMap<Integer, TreeSet<Regulator>>();
        for(int g = 1; g <= this.numGenes; g++) {
            this.regulators.put(g, new TreeSet<Regulator>());
        }
        
        this.corrs = new double[numGenes][numGenes];
        this.C = new int[numGenes][numGenes];
        this.fixedGenes = new boolean[numGenes];
                
        for (int rg = 0; rg < numGenes; rg++) {
            this.fixedGenes[rg] = false;
            
            for (int tg = 0; tg < numGenes; tg++) {
                if (rg != tg) {
                    C[rg][tg] = Regulator.RANK_UNKNOWN;
                } else {
                    C[rg][tg] = Regulator.RANK_NOLINK;
                }
                
                this.corrs[rg][tg] = 0;
            }
        }
    }
    
    public void setData(boolean[][] values) {
        this.data = values;
    }
    
    public void setExpGenes(int[] values) {
        this.expGenes = values;
    }
    
    public void setWildTypes(int[] values) {
        this.wildTypes = values;
    }
    
    public boolean[][] getData() {
        return this.data;
    }
    
    public int[] getExpGenes() {
        return this.expGenes;
    }
    
    public int[] getWildTypes() {
        return this.wildTypes;
    }
    
    public HashMap<Integer, TreeSet<Regulator>> getRegulators() {
        return this.regulators;
    }        
    
    public void findRegulators() {
        boolean[][] colValues = new boolean[this.numGenes][];
        
        for (int g = 1; g <= this.numGenes; g++) {
            colValues[g - 1] = Utils.parse(data, g - 1);
        }            
        double[][] d_colValues = Utils.toDouble(colValues);
        
        //correlation
        for (int g = 1; g <= this.numGenes; g++) {
            for (int eg = g; eg <= this.numGenes; eg++) {
                if (g != eg) {
                    double[] results = Utils.calCorrelation(d_colValues[g - 1], d_colValues[eg - 1],
                            d_colValues[g - 1].length, 2);
                    
                    if(Double.isNaN(results[0])) {
                        results[0] = 0;
                    }
                    this.corrs[eg - 1][g - 1] = results[0];
                    this.corrs[g - 1][eg - 1] = results[0];
                } else {
                    this.corrs[eg - 1][g - 1] = 1;
                }
            }
        }
        
        for(int g = 1; g <= this.numGenes; g++) {
            findRegulators(g);
        }
        
        //calculate Constraint matrix C
        for (int g = 1; g <= this.numGenes; g++) {
            TreeSet<Regulator> regGenes = regulators.get(g);
            
            for (Regulator reg: regGenes) {
                this.C[reg.regulator - 1][g - 1] = reg.rank;
            }
        }
    }
    
    private void findRegulators(int eg) {
        //values[0]: first row contains wild-type values of all genes
        //values[mutatedGene]: this row contains gene expression values of all genes 
        //                      as the 'mutatedGene' is perturbed
        int ko;
        TreeSet<Regulator> links = new TreeSet<Regulator>();
        int numExps = data.length;
        
        System.out.printf("<Gene %d>\t regulators = ", eg);
        for(int ex = 0; ex < numExps; ex++) {
            int pertG = expGenes[ex];
            
            if(pertG < 0) continue;
            
            if(pertG != eg) {
                int wtIndex = wildTypes[ex];
                
                if(InferBN.MODE_STRICT) {
                    if(data[ex][pertG - 1] == data[wtIndex][pertG - 1]) continue;
                }
                
                if(data[ex][pertG - 1] == false) {
                    //Knockout case
                    ko = -1;
                } else {
                    ko = 1;
                }
                                
                int egChange = Utils.minus(data[ex][eg - 1], data[wtIndex][eg - 1]);
                int type = -1;
                if(egChange != 0) {
                    if(ko * egChange > 0) {                        
                        type = 1;
                    } 
                    
                    double r = this.corrs[eg - 1][pertG - 1];
                    boolean sameSign = r * type > 0 ? true: false;
                    if(! InferBN.TF_TG_CORR_DIFF_SAMESIGN) sameSign = true;
                    
                    if(sameSign && Math.abs(r) >= InferBN.GENE_CORREL_THRESH_MID) {
                        if(type > 0) {
                            links.add(new Regulator(pertG, type, Regulator.RANK_DIRECT_POSITIVE));
                        } else {
                            links.add(new Regulator(pertG, type, Regulator.RANK_DIRECT_NEGATIVE));
                        }
                        
                        if(Math.abs(r) >= InferBN.GENE_CORREL_THRESH_MAX) {
                            this.fixedGenes[eg - 1] = true;
                            System.out.printf("<fix>");
                        }
                        System.out.printf("<dir>");
                        
                    } else {
                        if(InferBN.DATABASE != InferBN.ECOLI) {
                            if(r * type > 0 && Math.abs(r) >= InferBN.GENE_CORREL_THRESH_MIN) {
                                if(type > 0) {
                                    links.add(new Regulator(pertG, type, Regulator.RANK_INDIRECT_POSITIVE));
                                } else {
                                    links.add(new Regulator(pertG, type, Regulator.RANK_INDIRECT_NEGATIVE));
                                }
                            }
                            System.out.printf("<indir>");
                        }
                    }
                    
                    System.out.printf(pertG + " ");
                } else {
                    /*double r = this.corrs[eg - 1][g - 1];
                    if(Math.abs(r) < InferBN.GENE_CORREL_THRESH_MIN) {
                        links.add(new Regulator(g, 0, Regulator.RANK_NOLINK));
                        System.out.printf(g + "<No> ");
                    } else {
                        if(Math.abs(r) >= InferBN.GENE_CORREL_THRESH_MID) {
                            links.add(new Regulator(g, 0, Regulator.RANK_INDIRECT_NEUTRAL));
                            System.out.printf(g + "<Neu> ");
                        }
                    }*/
                }
            }
        }
        System.out.println();                                
        
        this.regulators.get(eg).addAll(links);
    }
    
    public void loadGoldData(String[][] data) {
        int numLines = data[0].length;
        
        for(int l = 0; l < numLines; l ++) {
            int srcGene = Integer.parseInt(data[0][l].substring(1));
            int tarGene = Integer.parseInt(data[1][l].substring(1));
            String strtype = data[2][l];
            int type = 1;
            if(strtype.equals("-")) type = -1;
            
            this.regulators.get(tarGene).add(new Regulator(srcGene, type, 0));
        }
    }
    
    public void outputRegulators(PrintWriter pw, String del) throws FileNotFoundException {
        //PrintWriter pw = new PrintWriter(new FileOutputStream(path),true);                
        
        pw.println("Target gene" + del + "Regulators");
        for(int g = 1; g <= this.numGenes; g++) {
            TreeSet<Regulator> dsGenes = this.regulators.get(g);
            pw.printf("%d", g);
            
            StringBuilder sbGenes = new StringBuilder();
            StringBuilder sbTypes = new StringBuilder();
            StringBuilder sbRanks = new StringBuilder();
            
            sbRanks.append("     Ranking");
            for(Regulator dg: dsGenes) {                                
                String type = "+";
                if(dg.interaction == -1) type = "-";
                sbTypes.append(del + type);
                
                sbGenes.append(del + dg.regulator + type);
                sbRanks.append(del + dg.rank);
            }                        
            
            pw.println(sbGenes.toString());
            //pw.println(sbTypes.toString());
            //pw.println(sbRanks.toString());
        }
                                             
        //pw.close();
    }
             
    public static GeneBData convert(ArrayList<Node> nodes) {
        GeneBData gdata = new GeneBData(nodes.size());

        for (Node node : nodes) {
            int targetID = Integer.valueOf(node.NodeID) + 1;
            int pos, Im, Om;
            LogicTable logic = node.getLogicTable();
            int numInputs = logic.input.size();

            for (int k = 0; k < numInputs; k++) {
                pos = logic.input.get(k);
                Im = logic.I.get(k);
                Om = logic.O.get(k);

                int type = 1;
                if (Im != Om) {
                    type = -1;
                }

                int srcID = Integer.valueOf(nodes.get(pos).NodeID) + 1;
                gdata.regulators.get(targetID).add(new Regulator(srcID, type, 0));
            }
        }
        return gdata;
    }
}