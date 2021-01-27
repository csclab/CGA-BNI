/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package bni.comp;

import bni.InferBN;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.TreeSet;

/**
 *
 * @author colin_PC2
 */
public class StatData {
    GeneBData gOrigin;
    GeneBData gdata;
    int TP = 0; //True Positve
    int FP = 0; //False Positve
    int TN = 0; //True Negative
    int FN = 0; //False Negative
    
    String sTP = ""; //True Positve
    String sFP = ""; //False Positve
    String sTN = ""; //True Negative
    String sFN = ""; //False Negative
    
    double Precision;
    double TPR;
    double Accuracy;        
    
    public StatData(GeneBData gOrigin, GeneBData gdata) {
        this.gOrigin = gOrigin;
        this.gdata = gdata;
    }
    
    public void stat() {
        int numGenes = this.gOrigin.numGenes;
        HashMap<Integer, TreeSet<Regulator>> origin_dsGenes = this.gOrigin.getRegulators();
        HashMap<Integer, TreeSet<Regulator>> predicted_dsGenes = this.gdata.getRegulators();        
        int[] type = {1, -1};
        
        if(InferBN.STAT_IGNORE_SIGN) {
            type = new int[]{0};
            System.out.println("[stat] Type length = " + type.length);
        }
        
        for(int tar_g = 1; tar_g <= numGenes; tar_g++) {
            TreeSet<Regulator> originLinks = origin_dsGenes.get(tar_g);
            TreeSet<Regulator> predictedLinks = predicted_dsGenes.get(tar_g);
            
            for(int reg_g = 1; reg_g <= numGenes; reg_g++) {
                if(reg_g == tar_g) continue;
                
                for(int t = 0; t < type.length; t++) {
                    Regulator link = new Regulator(reg_g, type[t], 0);
                    
                    if(predictedLinks.contains(link)) {
                        if (originLinks.contains(link)) {
                            ++TP;
                            sTP = sTP + link.getLink() + tar_g + ",";
                        } else {
                            ++FP;
                            sFP = sFP + link.getLink() + tar_g + ",";
                        }
                    } else {
                        if(originLinks.contains(link)) {
                            ++ FN;
                            sFN = sFN + link.getLink() + tar_g + ",";
                        } else {
                            ++ TN;
//                            sTN = sTN + link.getLink() + tar_g + ",";
                        }
                    }
                }
            }                                 
        }
        
        if(TP == 0) {
            this.Precision = 0;
            this.TPR = 0;
        } else {
            this.Precision = (double)TP / (TP + FP);
            this.TPR = (double)TP / (TP + FN);
        }
        
        this.Accuracy = (double)(TP + TN) / (TP + FP + TN + FN);        
    }
    
    public void output(String path, String del) throws FileNotFoundException {
        PrintWriter pw = new PrintWriter(new FileOutputStream(path),true);                
                
        this.gdata.outputRegulators(pw, del);        
        pw.println();
        
        pw.println("True Positive" + del + TP);
        pw.println("True Positive cases" + del + sTP);
               
        pw.println("False Positive" + del + FP);
        pw.println("False Positive cases" + del + sFP);
        
        pw.println("True Negative" + del + TN);
        pw.println("True Negative cases" + del + sTN);
        
        pw.println("False Negative" + del + FN);
        pw.println("False Negative cases" + del + sFN);
        
        pw.println("Precision" + del + Precision);
        pw.println("True Positive rate (Recall)" + del + TPR);
        pw.println("Accuracy" + del + Accuracy);
        pw.close();
    }
    
    public static void outputAverage(ArrayList<StatData> stats, String filename, String del) throws FileNotFoundException {
        double Precision = 0;
        double TPR = 0;
        double Accuracy = 0;

        PrintWriter pw = new PrintWriter(new FileOutputStream(filename), true);
        pw.println("ID" + del + "Precision" + del + "Recall" + del + "Accuracy");
        
        for (int i = 0; i < stats.size(); i++) {
            StatData stat = stats.get(i);
            pw.println(i + del + stat.Precision + del + stat.TPR + del + stat.Accuracy);
                    
            Precision += stat.Precision;
            TPR += stat.TPR;
            Accuracy += stat.Accuracy;
        }

        Precision = Precision / stats.size();
        TPR = TPR / stats.size();
        Accuracy = Accuracy / stats.size();
                
        pw.println("avg" + del + Precision + del + TPR + del + Accuracy);

        pw.close();
    }
    
}