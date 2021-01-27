/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package bni;

import static bni.InferBN.DIR_BASE;
import static bni.InferBN.DREAM3_MAX_NO_EDGES;
import static bni.InferBN.DREAM3_SIZE;
import static bni.InferBN.DREAM3_SIZE_NO;
import static bni.InferBN.DREAM3_SPECIES;
import static bni.InferBN.infer_DREAM;
import bni.comp.GeneBData;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.InputStream;
import mod.jmut.core.Calc;

/**
 *
 * @author colin_PC2
 */
public class ARACNE {
    public static final int ARACNE = 0;
    public static final int BC3NET = 1;
    public static final int BTR = 2;
    public static final int GENIE3 = 3;
    
    public static final double[] ARACNE_THRESH = {
        0.05, 0.05, 0.05
    };
    public static final String ARACNE_DIR = DIR_BASE + "Inference\\OtherTools\\ARACNE\\";
    public static final String BC3NET_DIR = DIR_BASE + "Inference\\OtherTools\\BC3NET\\";
    public static final String BTR_DIR = DIR_BASE + "Inference\\OtherTools\\BTR\\";
    public static final String GENIE3_DIR = DIR_BASE + "Inference\\OtherTools\\GENIE3\\";
    
    public static final String[] DIRS = {
        ARACNE_DIR, BC3NET_DIR, BTR_DIR, GENIE3_DIR
    };
    
    public static void run_ARACNE_DREAM() throws IOException, InterruptedException {
        String path = ARACNE_DIR + "DREAM3_TRANSPOSED\\";
        
        if(InferBN.DATABASE == InferBN.ECOLI) {
            path = ARACNE_DIR + "ECOLI_TRANSPOSED\\";
            
        } else if(InferBN.DATABASE == InferBN.RBN) {
            path = ARACNE_DIR + "RBN_TRANSPOSED\\";
        }
        
        File folder = new File(path);
        File[] listOfFiles = folder.listFiles();
        double thresh = 0.05;
        
        for (int i = 0; i < listOfFiles.length; i++) {
            if (listOfFiles[i].isFile()) {
                System.out.printf("File %d: %s\n", i + 1, listOfFiles[i].getName());
                
                if(listOfFiles[i].getName().contains("Size10_")) {
                    thresh = ARACNE_THRESH[0];
                } else if(listOfFiles[i].getName().contains("Size50_")) {                    
                    thresh = ARACNE_THRESH[1];
                } else if(listOfFiles[i].getName().contains("Size100_")) {                    
                    thresh = ARACNE_THRESH[2];                    
                } else if(listOfFiles[i].getName().contains("Ecoli_large")) {                    
                    thresh = 0.25;
                }
                
                if(InferBN.DATABASE == InferBN.RBN) {
                    thresh = 0.15;
                }
                
                String cmd = "java -jar " + ARACNE_DIR + "aracne2.jar -i " +
                            listOfFiles[i].getPath() + " -t " + thresh;
                Process proc = Runtime.getRuntime().exec(cmd);
                proc.waitFor();
                                
                InputStream in = proc.getInputStream();
//                InputStream err = proc.getErrorStream();
                byte b[] = new byte[in.available()];
                in.read(b, 0, b.length);
                System.out.println(new String(b));
            }
        }
    }
    
    public static void main(String[] args) throws FileNotFoundException, IOException, InterruptedException {
        
//        run_ARACNE_DREAM();    
        infer();
    }
    
    public static void infer() throws FileNotFoundException, IOException, InterruptedException {        
        
        int method = BC3NET;
        String dir_data = DIRS[method];
        
        if(InferBN.DATABASE == InferBN.DREAM3) {
            dir_data += "DREAM3_res\\";
        } else 
        if(InferBN.DATABASE == InferBN.GNW) {
            dir_data += "GNW_res\\";
        } else 
        if(InferBN.DATABASE == InferBN.ECOLI) {
            dir_data += "ECOLI_res\\";
        } else 
        if(InferBN.DATABASE == InferBN.RBN) {
            dir_data += "RBN_res\\";
        }
        
        for(int run = 1; run <= 1; ++run) {
            File theDir = new File("run" + run + "\\");
            
            if (!theDir.exists()) {
                try {
                    theDir.mkdir();
                    System.out.println("Created directory: " + theDir.getName());
                } catch (SecurityException se) {                    
                }
            }
        
            File folder = new File(dir_data + theDir.getName() + "\\");
            File[] listOfFiles = folder.listFiles();        

            for (int i = 0; i < listOfFiles.length; i++) {
                if (listOfFiles[i].isFile()) {
                    System.out.printf("File %d: %s\n", i + 1, listOfFiles[i].getName());

                    String[] params = listOfFiles[i].getName().split("_|\\.", -1);
                    int actSize = -1;
                    int sp = -1;
                    int numGenes = 0;
                    int maxNoEdges = Integer.MAX_VALUE;
                    String prexName = null;
                    String koFile = null;
                    String netFile = null;
                    
                    if(InferBN.DATABASE == InferBN.DREAM3) {
                        actSize = Utils.index(params[0], DREAM3_SIZE);
                        sp = Utils.index(params[1], DREAM3_SPECIES);
                        numGenes = DREAM3_SIZE_NO[actSize];
                        maxNoEdges = DREAM3_MAX_NO_EDGES[actSize];
                        
                    } else if(InferBN.DATABASE == InferBN.GNW) {
                        actSize = Utils.index(params[0], InferBN.GNW_SIZE);
                        
                        numGenes = Integer.parseInt(params[0].substring(4));
                        maxNoEdges = InferBN.GNW_MAX_NO_EDGES[actSize];
                        
                        prexName = params[0] + "_" + params[1];
                        String path = InferBN.DIR_GNW + params[0] + "_" + params[1].substring(0, 5) + "\\";
                        koFile = path + prexName + "_knockouts.tsv.bool.csv";
                        netFile = path + prexName + "_goldstandard_signed.tsv";
                        
                    } else if(InferBN.DATABASE == InferBN.ECOLI) {
                        numGenes = InferBN.E_COLI_NO_GENES;
                        maxNoEdges = InferBN.E_COLI_MAX_NO_EDGES;
                        
                        prexName = params[0] + "_" + params[1];
                        koFile = InferBN.E_COLI_DIR + InferBN.E_COLI_FILE;
                        netFile = InferBN.E_COLI_NET;
                        
                    } else if(InferBN.DATABASE == InferBN.RBN) {
                        
                        numGenes = Integer.parseInt(params[0].substring(4));
                        maxNoEdges = numGenes * 2;
                        
                        prexName = params[0] + "_" + params[1];
                        String path = InferBN.DIR_RBN + params[0] + "\\";
                        koFile = path + prexName + "_knockouts.tsv.bool.csv";
                        netFile = path + prexName + "_goldstandard_signed.tsv";
                        
                    }
                    
                    System.out.printf("Paras[] = %d, %d, %d, %d, %s, %s, %s\n", 
                            actSize, sp, numGenes, maxNoEdges, prexName, koFile, netFile);

                    GeneBData g_data = null;
                    if(method == ARACNE) {
                        g_data = Utils.load_ARACNE_results(listOfFiles[i].getPath(), numGenes, maxNoEdges);
                    } else 
                    if(method == BC3NET) {
                        g_data = Utils.load_BC3NET_results(listOfFiles[i].getPath(), numGenes);
                    } else 
                    if(method == BTR) {
                        g_data = bni.BTR.loadPredictedNetwork(listOfFiles[i].getPath(), numGenes);
                    } else 
                    if(method == GENIE3) {
                        g_data = Utils.load_GENEI3_results(listOfFiles[i].getPath(), 
                                       numGenes, maxNoEdges);
                    }

                    Calc.datas.clear();
                    System.gc();
                    infer_DREAM(theDir.getName(), actSize, sp, g_data, 
                            koFile, netFile, prexName, numGenes);
                }
            }//end for: files
        }
        
    }
}
