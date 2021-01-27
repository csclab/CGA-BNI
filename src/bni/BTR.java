/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package bni;

import static bni.InferBN.DREAM3_SIZE;
import static bni.InferBN.DREAM3_SIZE_NO;
import static bni.InferBN.DREAM3_SPECIES;
import static bni.InferBN.FILE_NET;
import bni.comp.GeneBData;
import bni.comp.Regulator;
import bni.comp.StatData;
import java.io.FileNotFoundException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.TreeSet;

/**
 *
 * @author colin_PC2
 */
public class BTR {
    public static final String BTR_DIR = "D:\\HCStore\\BioInformatics\\Others\\Inference\\OtherTools\\BTR\\DREAM3\\__DREAM_SIZE__\\";
    
    public static GeneBData loadPredictedNetwork(String path, int numGenes) {
        String [][] data = Utils.loadTextFile(path, ",", true);
        GeneBData gdata = new GeneBData(numGenes);
        
        int numLines = data[0].length;
        HashMap<Integer, TreeSet<Regulator>> regulators = gdata.getRegulators();
        
        for(int l = 0; l < numLines; l ++) {
            int tarGene = Integer.parseInt(data[2][l].substring(1));                        
            
            String strtype = data[1][l];
            int type = 1;
            if(strtype.equals("inhibits")) type = -1;
            
            String[] srcs = data[0][l].split("&");
            for(int g = 0; g < srcs.length; g++) {
                int srcGene = Integer.parseInt(srcs[g].substring(1));
            
                regulators.get(tarGene).add(new Regulator(srcGene, type, 0));
            }
        }
        
        gdata.finalSourceDest = true;
        return gdata;
    }
    
    public static void BTR_DREAM(int actSize) throws FileNotFoundException {        
        
        String outDir = BTR_DIR.replace("__DREAM_SIZE__", DREAM3_SIZE[actSize]);
        ArrayList<StatData> stats = new ArrayList<StatData>();
        
        for (int sp = 0; sp < DREAM3_SPECIES.length; sp++) {
            String pathOriginNetw = FILE_NET.replace("__DREAM_SIZE__", DREAM3_SIZE[actSize]);
            pathOriginNetw = pathOriginNetw.replace("__DREAM_SPECIES__", DREAM3_SPECIES[sp]);

            GeneBData gOrigin = Utils.loadOriginNetwork(pathOriginNetw, DREAM3_SIZE_NO[actSize]);
        
            String outPath = outDir + DREAM3_SPECIES[sp] + "_edges.csv";
            GeneBData gdata = loadPredictedNetwork(outPath, DREAM3_SIZE_NO[actSize]);
            
            StatData stat = new StatData(gOrigin, gdata);
            stat.stat();
            stat.output(outDir + DREAM3_SPECIES[sp] + "_metrics.csv", ",");
            
            stats.add(stat);
        }        
        
        StatData.outputAverage(stats, outDir + DREAM3_SIZE[actSize] + "_metrics_avg.csv", ",");
    }
    
    public static void main(String[] args) throws FileNotFoundException {
        int actSize = 0;
        
        BTR_DREAM(actSize);
    }
}
