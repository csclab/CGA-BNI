package bni;

import bni.comp.EdgeAddition2;
import bni.comp.GeneBData;
import bni.comp.NetInfo;
import bni.comp.Regulator;
import bni.comp.StatData;
import java.io.File;
import java.io.FileNotFoundException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedList;
import java.util.Queue;
import java.util.Random;
import java.util.Set;
import java.util.TreeSet;
import java.util.concurrent.CountDownLatch;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import mod.EdgeRemoval;
import mod.Knockout;
import mod.Mutation;
import mod.OverExpression;
import mod.jmut.core.Calc;
import mod.jmut.core.Config;
import mod.jmut.core.Util;
import mod.jmut.core.comp.Attractor;
import mod.jmut.core.comp.Interaction;
import mod.jmut.core.comp.LogicTable;
import mod.jmut.core.comp.NetData;
import mod.jmut.core.comp.Node;
import mod.jmut.core.comp.NodeInteraction;
import org.apache.commons.math3.distribution.EnumeratedIntegerDistribution;
import org.apache.commons.math3.random.JDKRandomGenerator;
import org.apache.commons.math3.random.RandomGenerator;

/**
 *
 * @author colin
 */
public class InferBN {        
    public static final String DELIM = ","; 
    public static final int BA_EDGE_TO_ADD = 1;
    public static final double BA_PROBABILITY = 0.5;
    
    public static boolean MODE_STRICT = false;//DREAM3//    true;//Ecoli//            false;//RBN//    
    public static int OFFSET_COL_EXP_DATA = 1;//DREAM3//  2;//Ecoli//     
    public static final boolean STAT_IGNORE_SIGN = true;
    public static boolean TF_TG_CORR_DIFF_SAMESIGN = !true;//RBN//  !true;//Ecoli//  true;//DREAM3//  
//    public static final boolean TF_TG_FIX_INDIRECT = true;
    
    public static final double GENE_CORREL_THRESH_MIN = 0.1;
    public static double GENE_CORREL_THRESH_MID = 0.25;//RBN//    0.05;//Ecoli//    0.5;//DREAM3//    
    public static final double GENE_CORREL_THRESH_MAX = 0.89;
    
    public static final double SCORE_PRECISION = Math.pow(10, -3);
//    public static final double GENE_SCORE_PRECISION = Math.pow(10, -2);
    public static final int SCORE_MAX_CONVERGE = 250;
    public static final int GENE_SCORE_NO_CONVERGE = 1000;
    
    public static final int SCORE_NO_NEIGHBORS = 1;
    public static final int SCORE_NO_TRIALS = 100;
    public static final int SCORE_NO_INITIAL_STATES = 2000;
    
    public static final int GA_MAX_ITERATION = 2000;
    public static final int GA_POPULATION_SIZE = 20;
    public static final double GA_ELITE_RATIO = 0.2;
    
    public static final String DIR_OUTPUT = "out\\";
    public static final String DIR_OUTPUT_RULE = DIR_OUTPUT + "rules\\";
    public static final String DIR_OUTPUT_STRUCTS = DIR_OUTPUT + "structs\\";
    
    public static final int DREAM3 = 0;
    public static final int GNW = 1;
    public static final int ECOLI = 2;
    public static final int RBN = 3;
    public static int DATABASE = RBN;
    
    public static final int[] DREAM3_SIZE_NO = {
        10, 50, 100
    };
    
    public static final String[] DREAM3_SIZE = {
        "Size10", "Size50", "Size100"
    };
    public static final String[] DREAM3_SPECIES = {
        "Ecoli1", "Ecoli2", "Yeast1", "Yeast2", "Yeast3"
    };
    
    public static final int[] DREAM3_MAX_NO_EDGES = {
        40, 200, 600
    };
    
    public static String DIR_BASE = "D:\\HCStore\\Papers\\InferGA_Manuscript\\data\\";
    public static String DIR_DATA = DIR_BASE + "Inference\\DREAM3 in silico challenge\\__DREAM_SIZE__\\DREAM3 data\\";
    public static String FILE_KO = "InSilico__DREAM_SIZE__-__DREAM_SPECIES__-null-mutants.tsv.bool.csv";
    public static String FILE_NET = DIR_BASE + "Inference\\DREAM3 in silico challenge\\__DREAM_SIZE__\\Networks\\InSilico__DREAM_SIZE__-__DREAM_SPECIES__.tsv";    
        
    //--- E. coli section (M3D & regulonDB)
    public static String E_COLI_DIR = DIR_BASE + "Inference\\E_coli_v4_Build_6\\";
    public static String E_COLI_FILE = "E_coli_v4_Build_6_exps466_ko.csv.bool.csv";
    public static String E_COLI_NET = E_COLI_DIR + "network\\E_coli_M3D_RegulonDB.txt";
    public static final int E_COLI_NO_GENES = 925;  //947;  //1424;
    public static final int E_COLI_MAX_NO_EDGES = 3200;
    
    //--- GNW section
    public static final String DIR_GNW = DIR_BASE + "Inference\\GNW_data\\";
    public static final String[] GNW_SIZE = {
        "Size10", "Size50", "Size100", "Size250", "Size500"
    };
    public static final int[] GNW_MAX_NO_EDGES = {
        30, 140, 300, 750, 1150
    };
    
    //--- RBN section
    public static String DIR_RBN = DIR_BASE + "Inference\\RBN_data\\";
    
    
    public static void main(String[] args) throws FileNotFoundException, InterruptedException {
        
        for(int i=0; i < args.length; i++) {
        
            System.out.println(args[i]);
        }
        
        DIR_BASE = args[0];
        DIR_DATA = DIR_BASE + "Inference\\DREAM3 in silico challenge\\__DREAM_SIZE__\\DREAM3 data\\";
        FILE_KO = "InSilico__DREAM_SIZE__-__DREAM_SPECIES__-null-mutants.tsv.bool.csv";
        FILE_NET = DIR_BASE + "Inference\\DREAM3 in silico challenge\\__DREAM_SIZE__\\Networks\\InSilico__DREAM_SIZE__-__DREAM_SPECIES__.tsv";    

        //--- E. coli section (M3D & regulonDB)
        E_COLI_DIR = DIR_BASE + "Inference\\E_coli_v4_Build_6\\";
        E_COLI_FILE = "E_coli_v4_Build_6_exps466_ko.csv.bool.csv";
        E_COLI_NET = E_COLI_DIR + "network\\E_coli_M3D_RegulonDB.txt";
        
        DIR_RBN = DIR_BASE + "Inference\\RBN_data\\";
    
        String typeExperiment = args[1];
        int noTrials = Integer.valueOf(args[2]);
        
        if(typeExperiment.equalsIgnoreCase("ecoli")) {
            MODE_STRICT = true;
            OFFSET_COL_EXP_DATA = 2;
            GENE_CORREL_THRESH_MID = 0.05;
            DATABASE = ECOLI;
        }
        
        if(typeExperiment.equalsIgnoreCase("dream")) {
            TF_TG_CORR_DIFF_SAMESIGN = true;
            GENE_CORREL_THRESH_MID = 0.5;
            DATABASE = DREAM3;
        }
        
        if(typeExperiment.equalsIgnoreCase("rbn")) {
            OFFSET_COL_EXP_DATA = 2;
            DATABASE = RBN;
        }
        
        
        execute(noTrials);
    }
    
    public static void execute(int noTrials) throws FileNotFoundException, InterruptedException {
        
        for(int run = 1; run <= noTrials; ++run) {
            File theDir = new File("run" + run + "\\");
            
            if (!theDir.exists()) {
                try {
                    theDir.mkdir();
                    System.out.println("Created directory: " + theDir.getName());
                } catch (SecurityException se) {                    
                }
            }
            
            if(DATABASE == DREAM3) {
            for (int actSize = 0; actSize < DREAM3_SIZE.length; actSize++) {
                for (int sp = 0; sp < DREAM3_SPECIES.length; sp++) {
                    Calc.datas.clear();
                    System.gc();
                    GeneBData g_data = null;
                    
//                    g_data = Utils.load_GENEI3_results(DIR_DATA_GENIE3 + theDir.getName() + "\\"
//                                       + DREAM3_SIZE[actSize] + "_" + DREAM3_SPECIES[sp] + ".txt", 
//                                       DREAM3_SIZE_NO[actSize], DREAM3_MAX_NO_EDGES[actSize]);

                    infer_DREAM(theDir.getName(), actSize, sp, g_data, null, null, null, -1);
                }
            }
            } else if(DATABASE == GNW || DATABASE == RBN) {
                //GNW or RBN data
                String dir = DIR_GNW;
                if(DATABASE == RBN) dir = DIR_RBN;
                
                File folder = new File(dir);
                File[] subFolders = folder.listFiles();                
        
                for (int i = 0; i < subFolders.length; i++) {
                    if (! subFolders[i].isDirectory()) continue;

                    File[] files = subFolders[i].listFiles();

                    for (int f = 0; f < files.length; f++) {
                        if (! files[f].isFile()) continue;                                                   

                        if(! files[f].getName().contains("_knockouts.tsv.bool.csv")) continue;
                        System.out.printf("KO File %d: %s\n", f + 1, files[f].getName());
                        
                        Calc.datas.clear();
                        System.gc();
                        GeneBData g_data = null;

                        String filename_net = files[f].getName().replace("_knockouts.tsv.bool.csv", 
                                "_goldstandard_signed.tsv");
                        String[] params = filename_net.split("_", -1);
                        int numGenes = Integer.parseInt(params[0].substring(4));
                        String prexName = params[0] + "_" + params[1];
                        
                        infer_DREAM(theDir.getName(), -1, -1, g_data, files[f].getPath(),
                                subFolders[i].getPath() + "\\" + filename_net, prexName, numGenes);
                    }
                }
            } else if(DATABASE == ECOLI) {
                Calc.datas.clear();
                System.gc();
                GeneBData g_data = null;

                int numGenes = E_COLI_NO_GENES;
                String prexName = "Size1424_Ecoli";

                infer_DREAM(theDir.getName(), -1, -1, g_data, E_COLI_DIR + E_COLI_FILE,
                        E_COLI_NET, prexName, numGenes);
            }
        }
//        infer_ECOLI();
    }
    
    public InferBN(GeneBData g_initial_data) {
        
        rand_pert_type = new EnumeratedIntegerDistribution(this.perts_ToGenerate, this.perts_Prob);
                
        rand_select_tarGene = new Random();
        rand_select_input = new JDKRandomGenerator();
        rand_select_member = new JDKRandomGenerator();
        
        this.g_initial_data = g_initial_data;
        ExaminingRules = new HashSet<String>();
        
        this.fixedInputs = new int[this.g_initial_data.numGenes][];
        HashMap<Integer, TreeSet<Regulator>> regulators = g_initial_data.getRegulators();
        
        for (int g = 1; g <= this.g_initial_data.numGenes; g++) {
            TreeSet<Regulator> regGenes = regulators.get(g);
            int noFixedInputs = Utils.countFixedRegulators(regGenes);
            this.fixedInputs[g - 1] = new int[noFixedInputs];
            
            int rg = 0;
            for (Regulator reg: regGenes) {
                if(reg.rank == Regulator.RANK_DIRECT_POSITIVE
                        || reg.rank == Regulator.RANK_DIRECT_NEGATIVE) {
                    this.fixedInputs[g - 1][rg ++] = reg.regulator - 1;
                }
            }
        }
        Utils.printFixedInputs(this.fixedInputs);
        
//        int numExperiments = g_initial_data.getData().length;
        this.gScore_thresh = 0;//(int)Math.floor(GENE_SCORE_PRECISION * g_initial_data.numGenes * numExperiments);
//        Config.out("init", "gScore_thresh = " + gScore_thresh);
        
        this.prev_g_scores = new double[this.g_initial_data.numGenes];
        this.g_noConverges = new int[this.g_initial_data.numGenes];
        this.fixedGenes = this.g_initial_data.fixedGenes;
        
        for (int g = 0; g < this.g_initial_data.numGenes; g++) {
            this.prev_g_scores[g] = -1;
            this.g_noConverges[g] = 0;
//            this.fixedGenes[g] = false;
        }
    }
    
    public static void infer_DREAM(String workDir, int actSize, int sp, 
            GeneBData g_initial_data, 
            String fileKO, String pathOriginNetw, String prexName, int numGenes) 
            throws FileNotFoundException, InterruptedException {                        
        
        if(DATABASE == DREAM3) {           
            String path = DIR_DATA.replace("__DREAM_SIZE__", DREAM3_SIZE[actSize]);
            fileKO = FILE_KO.replace("__DREAM_SIZE__", DREAM3_SIZE[actSize]);
            fileKO = fileKO.replace("__DREAM_SPECIES__", DREAM3_SPECIES[sp]);
            fileKO = path + fileKO;

            pathOriginNetw = FILE_NET.replace("__DREAM_SIZE__", DREAM3_SIZE[actSize]);
            pathOriginNetw = pathOriginNetw.replace("__DREAM_SPECIES__", DREAM3_SPECIES[sp]);

            numGenes = DREAM3_SIZE_NO[actSize];
            prexName = DREAM3_SIZE[actSize] + "_" + DREAM3_SPECIES[sp];
        }
        
        String[] dirs = createFolders(workDir, prexName + "_" + System.currentTimeMillis());
        GeneBData gOrigin = Utils.loadOriginNetwork(pathOriginNetw, numGenes);
        //gOrigin.outputDownStream(paras[0] + "origin.csv", ",");
                    
        if(g_initial_data == null) {
            g_initial_data = get_initialData(fileKO, numGenes, OFFSET_COL_EXP_DATA);        
        } else {
            String[][] data = Utils.loadTextFile(fileKO, ",", false);
        
            g_initial_data.setData(Utils.parseBoolData(data, OFFSET_COL_EXP_DATA, data.length - 1));        
            g_initial_data.setExpGenes(Utils.parseIntData(data, 0));
            if(DATABASE == DREAM3) {
                g_initial_data.setWildTypes(Utils.makeIntArray(data[0].length, 0));
            } else {
                g_initial_data.setWildTypes(Utils.parseIntData(data, 1));        
            }
            
            Config.out("infer_DREAM", "Retrieved experiment data!");
        }
        
        Config.out("infer_DREAM", "Searching Boolean networks ...");
        infer(g_initial_data, gOrigin, dirs, prexName);
    }
    
    public static void infer_ECOLI(String workDir) throws FileNotFoundException, InterruptedException {
        int numGenes = E_COLI_NO_GENES;
        String[] dirs = createFolders(workDir, "Size1424_Ecoli_" + System.currentTimeMillis());                

        GeneBData gOrigin = Utils.loadOriginNetwork(E_COLI_NET, numGenes);
        //gOrigin.outputDownStream(paras[0] + "origin.csv", ",");

        GeneBData g_initial_data = get_initialData(E_COLI_DIR + E_COLI_FILE, numGenes, OFFSET_COL_EXP_DATA);        
        
        Config.out("infer_ECOLI", "Searching Boolean networks ...");
        infer(g_initial_data, gOrigin, dirs, "Size1424_Ecoli");
    }
    
    public static void infer(GeneBData g_initial_data, GeneBData gOrigin, 
            String[] dirs, String prexName) throws FileNotFoundException, InterruptedException {
        
        InferBN infBN = new InferBN(g_initial_data);
        long startTime = System.currentTimeMillis();
//        ArrayList<NetInfo> optimalNets = infBN.searchBN(1024, 
//                SCORE_MAX_ITERATION, SCORE_NO_CONVERGE, SCORE_NO_NEIGHBORS);
        ArrayList<NetInfo> optimalNets = infBN.GA_searchBN(SCORE_NO_INITIAL_STATES, 
                GA_MAX_ITERATION, SCORE_MAX_CONVERGE, GA_POPULATION_SIZE, GA_ELITE_RATIO);
        if(optimalNets.isEmpty()) return;
        
        long searchTime = (System.currentTimeMillis() - startTime) / 1000;
        Config.out("infer", "Finished searching in " + searchTime + " seconds.");
        Utils.outputDynamicsAccuracy(dirs[0] + prexName,
                optimalNets, searchTime, DELIM);
        
        ArrayList<StatData> stats = new ArrayList<StatData>();
        for(NetInfo net: optimalNets) {
            Node.createRules(net.netD.nodes);
            Utils.outputRules(dirs[1], net.netD);
            
            GeneBData gdata = GeneBData.convert(net.netD.nodes);
            StatData stat = new StatData(gOrigin, gdata);
            stat.stat();
            stat.output(dirs[2] + prexName 
                    + "_struct_" + net.netD.networkName + ".csv", ",");
            stats.add(stat);
        }                
        StatData.outputAverage(stats, dirs[0] + prexName 
                    + "_struct_avg.csv", ",");
        
    }
    
    public static GeneBData get_initialData(String filePath,
            int numGenes, int offCol) 
            throws FileNotFoundException {                                       
        
        System.out.println("Working directory: " + filePath);
        
        GeneBData gdata = new GeneBData(numGenes);
                
        String[][] data = Utils.loadTextFile(filePath, ",", false);
        
        gdata.setData(Utils.parseBoolData(data, offCol, data.length - 1));        
        gdata.setExpGenes(Utils.parseIntData(data, 0));        
        if (DATABASE == DREAM3) {
            gdata.setWildTypes(Utils.makeIntArray(data[0].length, 0));
        } else {
            gdata.setWildTypes(Utils.parseIntData(data, 1));
        }
                            
        gdata.findRegulators();
        
        System.out.println("Found initial regulators successfully!");
        return gdata;
    }                      
    
    public ArrayList<NetInfo> searchBN(int noStates, 
            int maxIteration, int noConverge, int noNeighbors) {
        int numGenes = g_initial_data.numGenes;
//        boolean[][] data = g_initial_data.getData();
        ArrayList<NetInfo> optimalNets = new ArrayList<NetInfo>();
        
        int noNet = 1;
        String netName = "RBN" + String.valueOf(noNet);
        NetInfo initNet = this.createRBN_BA(netName, numGenes, BA_EDGE_TO_ADD, 
                BA_PROBABILITY, g_initial_data);        
        Config.out("searchBN", "Created an initial random boolean network!");
                
        String stateSet = NetData.generateInitialStates(netName, String.valueOf(noStates));                        
        
        this.refineNetwork(initNet.netD);
        Utils.printInputs(initNet.netD.nodes);
        
        Queue<NetInfo> Q = new LinkedList<NetInfo>();
        Q.add(initNet);
        double minScore = Double.MAX_VALUE;        
        int noIter = 0;
        boolean optimal;
        
        while (!Q.isEmpty()) {
            ++ noIter;
            
            NetInfo net = Q.remove();
//            System.out.println("Processing the network: " + net.netD.networkName + " -------------------------------------------------------------------->>>");
            
//            boolean succ = this.calScore(net, stateSet, this.g_initial_data);
//            if (succ == false) {
//                Config.out("searchBN", "Failed to calculate the score of the network!");
//                continue;
//            }
            
            Config.out("searchBN", "Network " + net.netD.networkName + " - avg Score = " + net.avg_score);
            optimal = true;
            if(Math.abs(minScore - net.avg_score) <= SCORE_PRECISION) {
                optimalNets.add(net);
                if(optimalNets.size() >= noConverge) {
                    Config.out("searchBN", "The scores of optimal networks are converged! Stoped!");
                    break;
                }
            } else {
                if(minScore - net.avg_score > SCORE_PRECISION) {
                    minScore = net.avg_score;
                    optimalNets.clear();
                    optimalNets.add(net);
                } else {
                    //score of "net" is bad
                    optimal = false;
                }
            }
            
            if(net.avg_score == 0) {
                Config.out("searchBN", "Found a perfect network with ZERO score! Stoped!");
                break;
            }
            
            //Remove the network from Calc
            Calc.datas.remove(net.netD.networkName);
            if(noIter >= maxIteration) {
                Config.out("searchBN", "The maximum number of iterations is reached! Stoped!");
                break;                
            }                                    
            
            //Update number of convergences for each gene
            for (int g = 0; g < this.g_initial_data.numGenes; g++) {
                if(Double.compare(net.g_scores[g], this.prev_g_scores[g]) != 0) {
                    this.g_noConverges[g] = 0;
                } else {
                    this.g_noConverges[g] += 1;
                }
                
                this.prev_g_scores[g] = net.g_scores[g];
            }            
            
            ArrayList<NetInfo> neighborNets = new ArrayList<NetInfo>();
            if(optimal == true) {
                //Perturb the network to get newly modified networks
                neighborNets.addAll(this.getNeighborNetworks(net, noNeighbors, noNet + 1));                
            } else {
                //Stop perturb the network with a bad score
                if(Q.isEmpty()) {
                    //find neighbors from all optimal networks
                    int _noNet = noNet;
                    for (NetInfo opnet: optimalNets) {
                        ArrayList<NetInfo> _neighs = this.getNeighborNetworks(opnet, noNeighbors, _noNet + 1);
                        neighborNets.addAll(_neighs);
                        _noNet += _neighs.size();
                    }                    
                }
            }                                    
            
            for (int i = 0; i < neighborNets.size(); i++) {
                Q.add(neighborNets.get(i));
                Calc.datas.put(neighborNets.get(i).netD.networkName, neighborNets.get(i).netD);
            }
            noNet += neighborNets.size();
        }
        
        if(Q.isEmpty()) {
            Config.out("searchBN", "The queue is empty! Stoped!");
        }
        
        Config.out("searchBN", "Found " + optimalNets.size() + " optimal networks with score = " + minScore);
        return optimalNets;
    }
    
    public NetInfo createRBN_BA(String netName, int numGenes, int edgesToAdd, 
            double probability, GeneBData gdata) {
        ArrayList<Interaction> inatemp = new ArrayList<Interaction>();
        int[][] posPaths = new int[numGenes][numGenes];
        int[][] negPaths = new int[numGenes][numGenes];
        
        for (int rg = 0; rg < numGenes; rg++) {
            for (int tg = 0; tg < numGenes; tg++) {                                
                posPaths[rg][tg] = 0;
                negPaths[rg][tg] = 0;                
            }
        }
        
        try {
            RandomGenerator random = new JDKRandomGenerator();
            int degrees[] = new int[numGenes];

            int numofedges = 0;                        
            int inatype;
            int[][] C = gdata.C;
            
            if(gdata.finalSourceDest) {
                HashMap<Integer, TreeSet<Regulator>> regulators = gdata.getRegulators();
                
                for (int g = 1; g <= gdata.numGenes; g++) {
                    TreeSet<Regulator> regGenes = regulators.get(g);
                    
                    for (Regulator reg: regGenes) {
                        inatemp.add(numofedges, new Interaction(Integer.toString(reg.regulator - 1), 
                                    1, Integer.toString(g - 1)));                            
                        numofedges++;
                    }
                }
                System.out.println("[createRBN_BA] numofedges = " + numofedges);
                
            } else {//non-final source dest
                            
            //direct links
            for (int rg = 0; rg < numGenes; rg++) {
                for (int tg = 0; tg < numGenes; tg++) {
                    if (rg == tg) continue;
                    
                    switch(C[rg][tg]) {
                        case Regulator.RANK_DIRECT_POSITIVE:
                        case Regulator.RANK_DIRECT_NEGATIVE:
                            inatype = 1;
                            if (C[rg][tg] == Regulator.RANK_DIRECT_NEGATIVE) {
                                inatype = -1;
                            }

                            inatemp.add(numofedges, new Interaction(Integer.toString(rg), inatype, Integer.toString(tg)));
                            degrees[rg]++;
                            degrees[tg]++;
                            numofedges++;
                            
                            Utils.modifyLink(rg, inatype, tg, C, posPaths, negPaths, true, Regulator.MODE_ADD);
                            break;
                            
                        default:
                            break;
                    }
                }
            }
            System.out.println("<BA> No. of direct edges = " + numofedges);
            
            //indirect links
            boolean succ;
            for (int rg = 0; rg < numGenes; rg++) {
                for (int tg = 0; tg < numGenes; tg++) {
                    if (rg == tg) continue;
                    if(Utils.pairOK(rg, tg, C, posPaths, negPaths) == -1) continue;
                    
                    switch(C[rg][tg]) {
                        case Regulator.RANK_INDIRECT_POSITIVE:
                            succ = Utils.modifyLink(rg, 1, tg, C, posPaths, negPaths, false, Regulator.MODE_ADD);
                            if(succ) {
                                inatemp.add(numofedges, new Interaction(Integer.toString(rg), 1, Integer.toString(tg)));
                                degrees[rg]++;
                                degrees[tg]++;
                                numofedges++;
                            }
                            break;
                        
                        case Regulator.RANK_INDIRECT_NEGATIVE:
                            succ = Utils.modifyLink(rg, -1, tg, C, posPaths, negPaths, false, Regulator.MODE_ADD);
                            if(succ) {
                                inatemp.add(numofedges, new Interaction(Integer.toString(rg), -1, Integer.toString(tg)));
                                degrees[rg]++;
                                degrees[tg]++;
                                numofedges++;
                            }
                            break;
                                
                        case Regulator.RANK_INDIRECT_NEUTRAL:
                            inatype = 1;
                            succ = Utils.modifyLink(rg, inatype, tg, C, posPaths, negPaths, false, Regulator.MODE_ADD);
                            if(succ == false) {
                                inatype = -1;
                                succ = Utils.modifyLink(rg, inatype, tg, C, posPaths, negPaths, false, Regulator.MODE_ADD);
                            }
                            if(succ) {
                                inatemp.add(numofedges, new Interaction(Integer.toString(rg), inatype, Integer.toString(tg)));
                                degrees[rg]++;
                                degrees[tg]++;
                                numofedges++;
                            }
                            
                            break;
                            
                        default:
                            break;
                    }
                }
            }
            System.out.println("<BA> No. of direct + indirect edges = " + numofedges);
            
            if(! (DATABASE == ECOLI && numofedges >= E_COLI_MAX_NO_EDGES)) {
            //add remaining links
            int[] index_select = new int[numGenes];
            double[] probs_select = new double[numGenes];
            boolean[] visit = new boolean[numGenes];
            
            for(int i = 0; i < numGenes; i++) {
                index_select[i] = i;
            }
            
            for (int tg = 0; tg < numGenes; tg++) {
                if(degrees[tg] > 0) continue;
                if(gdata.fixedGenes[tg] == true) continue;
                                
                double oldTotalDegrees = 2.0d * numofedges;
                int cnt = 0;
                
                for (int i = 0; i < numGenes; i++) {
                    probs_select[i] = degrees[i] / oldTotalDegrees;
                    visit[i] = false;
                    
                    if(probs_select[i] == 0 || oldTotalDegrees == 0) {
                        probs_select[i] = 0.5;
                    }
                }
                visit[tg] = true;
                EnumeratedIntegerDistribution rand = new EnumeratedIntegerDistribution(random, index_select, probs_select);
        
                while (cnt < edgesToAdd) {
                    int rg = rand.sample();
                    if (visit[rg] == true) {
                        continue;
                    }

                    inatype = (Math.random() < probability) ? -1 : 1;
                    Interaction temp = new Interaction(Integer.toString(rg), inatype, Integer.toString(tg));

                    if (Utils.checkExistInteraction(temp, numofedges, inatemp) == false) {
                        succ = Utils.modifyLink(rg, inatype, tg, C, posPaths, negPaths, false, Regulator.MODE_ADD);
                        if (succ == false) {
                            inatype = -inatype;
                            succ = Utils.modifyLink(rg, inatype, tg, C, posPaths, negPaths, false, Regulator.MODE_ADD);
                        }

                        if (succ == true) {
                            temp.InteractionType = inatype;
                            inatemp.add(numofedges, temp);
                            degrees[rg]++;
                            degrees[tg]++;
                            numofedges++;
                            
                            ++ cnt;
                        }
                    }

                    visit[rg] = true;
                    if(! Utils.exist(false, visit)) break;
                }
            }
            }
            
            }

        } catch (Exception e) {
            e.printStackTrace();
        }

//        for(Interaction ind: inatemp) {
//            System.out.println(ind.getID());
//        }
        
        NetData.loadNetwork_v2(netName, inatemp, numGenes);
        NetData.generateCurrentRule(netName, 0);
        NetInfo initNet = new NetInfo(Calc.datas.get(netName), posPaths, negPaths);
        
        return initNet;
    }    

    private static boolean[][] toBools(Attractor att, int numGenes) {        
        boolean[][] data = new boolean[att.Length][numGenes];
                
        for (int i = 0; i < att.Length; i++) {
            String state = att.States.get(i);
            
            for (int g = 0; g < numGenes; g++) {
                data[i][g] = state.charAt(g) == '1'? true: false;
            }
        }
        
        return data;
    }
    
    public static boolean calScore(NetInfo net, boolean[][] states, 
            boolean[][] base_data, int[] expGenes, int[] wildTypes) {        
//        long startTime = System.currentTimeMillis();
        
        Calc cal = new Calc();
//        boolean[][] base_data = g_initial_data.getData();
//        int[] expGenes = g_initial_data.getExpGenes();
        int numGenes = base_data[0].length;
        
        //Find all wild-type attractors
        ArrayList<Attractor> AllAttractors = cal.findAttractors_v2(net.netD, states, false);
        
//        System.out.printf("Time to find all wild-type attractors = %.2f\n", 
//                (double)(System.currentTimeMillis() - startTime) / 1000);
//        System.out.println("Number of wild-type attractors = " + Calc.AllAttractors.size());
        
        if(AllAttractors == null) {
            Config.out("calScore", "Failed to find all wild-type attractors of the network!");
            return false;
        }
        
//        for(Attractor att: Calc.AllAttractors) {
//            Util.printAttractor(att);
//        }
        
        //Check if exists an attractor same with the wild-type one in the true attractors (from expression data)
        int numExps = expGenes.length;
        ArrayList<Attractor> wtAttractors = new ArrayList<Attractor>();
                
        for(int ex = 0; ex < numExps; ex++) {
            int pertG = expGenes[ex];
            
            if(pertG != -1) break;
            
            double minScore = Double.MAX_VALUE;
            Attractor wt_att = null;

            for(int i = 0; i < AllAttractors.size(); i++) {
                Attractor att = AllAttractors.get(i);
//                if(att.Length > 1) continue;

                boolean[][] data = InferBN.toBools(att, numGenes);
                double[] g_scores = InferBN.compare(data, base_data[ex], numGenes);
                double score = Utils.sum(g_scores);
                if(score < minScore) {
                    minScore = score;
                    wt_att = att;
                }
            }
            
            if(wt_att != null) {
                wtAttractors.add(wt_att);
                
                //Calculate again the wild-type attractor score & add the score to the network
                boolean[][] data = InferBN.toBools(wt_att, numGenes);
                double[] __g_scores = InferBN.compare(data, base_data[ex], numGenes);
                
                if(ex == 0) {
                    net.g_scores = __g_scores;
                } else {
                    net.add(__g_scores);
                }
            }
        }
        
        if(wtAttractors.isEmpty()) {
            //cannot find the fix-point attractor => score = Double.MAX_VALUE
            Config.out("calScore", "Non-exist fix-point attractors => score = Double.MAX_VALUE!");
            return true;
        }
        
//        Config.out("calScore", "selected wild-type attractor: sum score = " + minScore);
//        Util.printAttractor(wt_att);        
        AllAttractors.clear();
                
        //Perturb the network like the expression data shown     
//        String wt_stateSet = NetData.setInitialState(net.netD.networkName, wt_att.States.get(0));
        ArrayList<boolean[][]> wtStates = new ArrayList<boolean[][]>();
        
        for(int i = 0; i < wtAttractors.size(); i++) {
            Attractor wt_att = wtAttractors.get(i);
            boolean[][] wt_state = Utils.convert_2DBool(wt_att.States.get(0));
            
            wtStates.add(wt_state);
        }
        
//        System.out.println("Finding the perturbed attractors ...");
//        startTime = System.currentTimeMillis();
        
        for (int ex = 0; ex < expGenes.length; ex++) {
            int pertG = expGenes[ex] - 1;
            
            if(pertG < 0) continue;            
            Mutation mut;
            int wtIndex = wildTypes[ex];
            boolean[][] wt_state = wtStates.get(wtIndex);
            
            //check knockout or overexpression mutation            
            if(base_data[ex][pertG] == true) {
                mut = new OverExpression();                
            } else {
                mut = new Knockout();
            }
            
            Util.perturbNode(net.netD.nodes.get(pertG), mut);
            ArrayList<Attractor> pertAttractors = cal.findAttractors_v2(net.netD, wt_state, false);
            if(pertAttractors == null) {
                Config.out("calScore", "Failed to find the perturbed attractor along with the mutation on gene " + pertG);
                return false;
            }            
            net.netD.nodes.get(pertG).restoreOriginalLogicTable();
            
            //Calculate the perturbed attractor score & add the score to the network
            Attractor pert_att = pertAttractors.get(0);
            boolean[][] data = InferBN.toBools(pert_att, numGenes);
            double[] g_scores = InferBN.compare(data, base_data[ex], numGenes);
            net.add(g_scores);
            
            pertAttractors.clear();
//            Config.out("calScore", "Perturbed gene " + pert_g + " : attractor with sum score = " + Utils.sum(g_scores));
//            Util.printAttractor(pert_att);
        }
        
        net.average(numGenes + 1);
        
//        System.out.printf("Time to find all perturbed attractors = %.2f\n", 
//                (double)(System.currentTimeMillis() - startTime) / 1000);
        
        return true;
    }
    
    private static double[] compare(boolean[][] data, boolean[] base_data, int numGenes) {
        double[] g_scores = new double[numGenes];
        java.util.Arrays.fill(g_scores, 0);

        for (int s = 0; s < data.length; s++) {
            for (int g = 0; g < numGenes; g++) {
                if (data[s][g] != base_data[g]) {
                    g_scores[g] += 1;
                }
            }
        }
        
        for (int g = 0; g < numGenes; g++) {
            g_scores[g] = g_scores[g] / data.length;
        }
        return g_scores;
    }
    
    public ArrayList<NetInfo> getNeighborNetworks(NetInfo net, int noNeighbors, int noNet) {
        ArrayList<NetInfo> neighborNets = new ArrayList<NetInfo>();        
        int cnt = 0;
        int noTrials = 0;
        //int noPerts = 1;
        double avg_score = Utils.sum(net.g_scores) / net.g_scores.length;
        int res;
        
        while(cnt < noNeighbors) {
            NetData newNetD = net.netD.createCopy();
            newNetD.networkName = "RBN" + String.valueOf(noNet + cnt);
            
            NetInfo newNet = new NetInfo(newNetD, net.posPaths, net.negPaths);
            //for(int p = 0; p < noPerts; p ++) {
                res = this.perturb(newNet, net.g_scores, avg_score);
                if(res == -1) break;
            //}
            
            //check if exists before
            String curRules = null;
            if(res == 1) {
                curRules = Node.getRulesInBinaryFormat(newNetD.nodes);
            }
            
            if(res == 0 || ExaminingRules.add(curRules) == false) {
                ++ noTrials;
//                Config.out("getNeighborNetworks", "Found a duplicate network. Retry again: " 
//                        /*+ noPerts*/ + "/" + noTrials);                
                
                if(noTrials > SCORE_NO_TRIALS) {
                    //++ noPerts;
                    //if(noPerts > net.g_scores.length) {
                    Config.out("getNeighborNetworks", "Failed to find a neighbor network.");    
                    break;
                    //} else {
                    //    noTrials = 0;
                    //    continue;
                    //}
                } else {
                    continue;
                }
            }
            
            //update nodes
            Utils.updateNodes(newNetD);
            
            neighborNets.add(newNet);
            noTrials = 0;
            //noPerts = 1;
            ++ cnt;
        }
        
        return neighborNets;
    }
    
    public int perturb(NetInfo net, double[] g_scores, double avg_score) {
        /*
         * return
         *  -1  cannot make a new network by mutation
         *  0   cannot make a new network by mutation, but can select another target gene
         *  1   OK
         */
        int numGenes = net.netD.nodes.size();
        
        //check if exists a gene scores > gScore_thresh
        boolean found = false;
        for(int g = 0; g < g_scores.length; g++){
            if(g_scores[g] > this.gScore_thresh 
                    && this.g_noConverges[g] < GENE_SCORE_NO_CONVERGE
                    && this.fixedGenes[g] == false
                    /*&& g_scores[g] >= avg_score*/) {
                found = true;
                break;
            }
        }        
        if(found == false) return -1;
        
        //randomly select a target gene
        int tarG;
        while(true) {
            tarG = rand_select_tarGene.nextInt(numGenes);
            if(g_scores[tarG] > this.gScore_thresh
                    && this.g_noConverges[tarG] < GENE_SCORE_NO_CONVERGE
                    && this.fixedGenes[tarG] == false
                    /*&& g_scores[tarG] >= avg_score*/) {
                break;
            }
        }
        
        /*
         * Select a type of perturbation:
         *      Change I                    (noInputs - noFixedInputs >= 1)
         *      Change O                    (noInputs - noFixedInputs >= 1)
         *      Change I & O                (noInputs >= 1)
         *      Swap two inputs ip1 & ip2   (noInputs - noFixedInputs >= 2)
         * 
         *      Remove one input            (noInputs - noFixedInputs >= 1)
         *      Insert one input            (noInputs < numGenes - 1)
         */
        LogicTable logic = net.netD.nodes.get(tarG).getLogicTable();
        int noInputs = logic.input.size();
        int noFixedInputs = this.fixedInputs[tarG].length;
        int pert_type;
        boolean done = false;
        
        double[][] corrs = this.g_initial_data.corrs;
        //edge-removal
        double[] probs_remove = new double[noInputs - noFixedInputs];
        int[] index_remove = new int[noInputs - noFixedInputs];
        /*double probs_sum = 0;
        
        for(int i = 0; i < noInputs - noFixedInputs; i++) {
            probs_sum += Math.abs(corrs[tarG][logic.input.get(i + noFixedInputs)]);
        }*/
        
        for(int i = 0; i < noInputs - noFixedInputs; i++) {
            probs_remove[i] = 1 - Math.abs(corrs[tarG][logic.input.get(i + noFixedInputs)]);// / probs_sum;
            index_remove[i] = i;
            
            if(probs_remove[i] == 0) {
                probs_remove[i] = 0.01;//fix bug: MathArithmeticException: array sums to zero
            }
        }
        
        //edge-addition
        ArrayList<Integer> temp = new ArrayList<Integer>();
        double[] probs_insert = null;
        int[] index_insert = null;
        boolean[] visit_insert = null;
        
        for(int i = 0; i < numGenes; i++) {
            if(i != tarG && Utils.exist(i, logic.input) == false
                    && Math.abs(corrs[tarG][i]) >= GENE_CORREL_THRESH_MIN) {
                temp.add(i);
            }
        }
        if(temp.size() > 0) {
            probs_insert = new double[temp.size()];
            index_insert = new int[temp.size()];
            visit_insert = new boolean[temp.size()];
            
            for(int i = 0; i < temp.size(); i++) {
                probs_insert[i] = Math.abs(corrs[tarG][temp.get(i)]);
                index_insert[i] = temp.get(i);
                visit_insert[i] = false;
            }
        }
        
        //check visit cases to avoid hang out
        if(noInputs < 1) {
            //noInputs >= 1: case 2 is always possible
            if(probs_insert == null) return 0;
        }
        
        while(! done) {
            pert_type = rand_pert_type.sample();                                
            int ip1, ip2, iSrc;                                
        
            switch(pert_type) {
                case 0:
                case 1:                
                    if(noInputs - noFixedInputs >= 1) {
                        //if(Utils.exist(logic.input.get(ip1), this.fixedInputs[tarG]) == false) break;
                        EnumeratedIntegerDistribution rand = new EnumeratedIntegerDistribution(rand_select_input, index_remove, probs_remove);
                        ip1 = noFixedInputs + rand.sample();

                        if(g_initial_data.finalSourceDest == false) {
                        int type = ( logic.I.get(ip1) == logic.O.get(ip1) ) ? 1: -1;
                        if(Math.abs(corrs[tarG][logic.input.get(ip1)]) >= GENE_CORREL_THRESH_MID) {
                            if(type * corrs[tarG][logic.input.get(ip1)] > 0) break;
                        }
                        
                        if(Utils.modifyLink(logic.input.get(ip1), type, tarG, 
                                g_initial_data.C, net.posPaths, net.negPaths, false, Regulator.MODE_CHG) == false) {
                            break;
                        }
                        }
                        
                        if(pert_type == 0) {
                            logic.I.set(ip1, 1 - logic.I.get(ip1));
                        } else {
                            logic.O.set(ip1, 1 - logic.O.get(ip1));
                            logic.O_default.set(1 - logic.O.get(noInputs - 1));
                        }
                        
                        //move the input "ip1" to the beginning of the non-fixed inputs series
                        logic.move(ip1, noFixedInputs);
                        done = true;
                    }
                    break;                         

                case 2:
                    if(noInputs >= 1) {
                        ip1 = rand_select_input.nextInt(noInputs);

                        logic.I.set(ip1, 1 - logic.I.get(ip1));
                        logic.O.set(ip1, 1 - logic.O.get(ip1));
                        logic.O_default.set(1 - logic.O.get(noInputs - 1));
                        
                        //move the input "ip1" to the beginning of the corresponding inputs series
                        if(ip1 < noFixedInputs) {
                            logic.move(ip1, 0);
                        } else {
                            logic.move(ip1, noFixedInputs);
                        }
                        done = true;
                    }
                    break;

                case 3:
                    if(noFixedInputs >= 2) {
                        for(int i = 0; i < noFixedInputs - 1; i++) {
                            ip1 = rand_select_input.nextInt(noFixedInputs);
                            ip2 = ip1;

                            while(ip2 == ip1) {
                                ip2 = rand_select_input.nextInt(noFixedInputs);
                            }
                            logic.swap(ip1, ip2);                        
                        }
                        done = true;
                    }
                    
                    if(noInputs - noFixedInputs >= 2) {
                        for(int i = 0; i < noInputs - noFixedInputs - 1; i++) {
                            ip1 = noFixedInputs + rand_select_input.nextInt(noInputs - noFixedInputs);
                            ip2 = ip1;

                            while(ip2 == ip1) {
                                ip2 = noFixedInputs + rand_select_input.nextInt(noInputs - noFixedInputs);
                            }
                            logic.swap(ip1, ip2);                        
                        }
                        done = true;
                    }
                    break;

                case 4:                    
                    //Remove one input            
                    if(g_initial_data.finalSourceDest == true) {
                        break;
                    }
                    
                    if(noInputs - noFixedInputs >= 1) {
                        EnumeratedIntegerDistribution rand = new EnumeratedIntegerDistribution(rand_select_input, index_remove, probs_remove);
                        ip1 = noFixedInputs + rand.sample();

                        iSrc = logic.input.get(ip1);
                        int type = ( logic.I.get(ip1) == logic.O.get(ip1) ) ? 1: -1;
                        
                        if(Utils.modifyLink(iSrc, type, tarG, 
                                g_initial_data.C, net.posPaths, net.negPaths, false, Regulator.MODE_DEL) == false) {
                            break;
                        }
                        
                        new EdgeRemoval().edgeMutation(iSrc, logic.input, logic.I, logic.O, logic.O_default);        
                        done = true;
                    }
                    break;

                case 5:
                    if(g_initial_data.finalSourceDest == true) {
                        break;
                    }
                    
                    //Insert one input, select a new input (also a source node)                
                    if(noInputs < numGenes - 1 && probs_insert != null) {                        
                        EnumeratedIntegerDistribution rand = new EnumeratedIntegerDistribution(
                                rand_select_input, index_insert, probs_insert);
                        iSrc = rand.sample();
                        
                        int index = Utils.index(iSrc, index_insert);
                        if(visit_insert[index] == true) break;
                        visit_insert[index] = true;
                                     
                        int type = 1;                        
                        if(Utils.modifyLink(iSrc, type, tarG, 
                                g_initial_data.C, net.posPaths, net.negPaths, false, Regulator.MODE_ADD) == false) {
                            type = -1;
                            if(Utils.modifyLink(iSrc, type, tarG, 
                                g_initial_data.C, net.posPaths, net.negPaths, false, Regulator.MODE_ADD) == false) {
                                break;
                            }
                        }
                        
                        Mutation mut = new EdgeAddition2(noFixedInputs, type);
                        mut.edgeMutation(iSrc, logic.input, logic.I, logic.O, logic.O_default);
                        done = true;
                    }
                    break;

                default:
                    break;
            }
            
            if(noInputs < 1 && ! Utils.exist(false, visit_insert)) {
                return 0;
            }
        }
        
        return 1;
    }
    
    private void refineNetwork(NetData netD) {
        for(int n = 0; n < netD.nodes.size(); n++) {
            LogicTable logic = netD.nodes.get(n).getLogicTable();
            int noInputs = logic.input.size();
            int[] fixed = this.fixedInputs[n];
            double[][] corrs = this.g_initial_data.corrs;
            
            if(fixed.length > 0) {
                for(int fi: fixed) {
                    if(! Utils.exist(fi, logic.input)) {
                        System.out.println("error in logic of node " + n + "_" + netD.nodes.get(n).NodeID);
                        System.out.println(java.util.Arrays.toString(fixed));
                        System.out.println(java.util.Arrays.toString(logic.input.toArray()));
                        
                        for(NodeInteraction ni: netD.in.get(netD.nodes.get(n).NodeID)) {
                            System.out.print(ni.Node + " ");
                        }
                        System.out.println();
                        break;
                    }
                }
            }
            
            if(fixed.length > 0) {
                //reorder fixed inputs to the head
                for(int i = 0; i < noInputs; i++) {
                    if(Utils.exist(logic.input.get(i), fixed) == false) {
                        boolean found = false;
                        
                        for(int j = i + 1; j < noInputs; j++) {
                            if(Utils.exist(logic.input.get(j), fixed) == true) {//is fixed input
                                logic.swap(i, j);
                                found = true;
                                break;
                            }
                        }
                        
                        if(found == false) {
                            break;
                        }
                    }
                }                
                logic.O_default.set(1 - logic.O.get(noInputs - 1));
                
                //remove all non-fixed inputs if exists a fixed input with high correlation
                /*boolean found = false;                
                
                for(int i = 0; i < fixed.length; i++) {
                    if(Math.abs(corrs[n][logic.input.get(i)]) >= GENE_CORREL_THRESH_MAX) {
                        found = true;
                        break;
                    }
                }*/
                
                if(this.fixedGenes[n]) {
                    for(int i = noInputs - 1; i >= fixed.length; i--) {
                        logic.remove(i);
                    }
                    logic.O_default.set(1 - logic.O.get(fixed.length - 1));
                    
                    //mark gene 'n' as a fixed target gene, and do not apply mutation
//                    this.fixedGenes[n] = true;
                }
            }
            
            //remove non-fixed inputs with very low correlation <= GENE_CORREL_THRESH_MIN
            noInputs = logic.input.size();
            for(int i = noInputs - 1; i >= fixed.length; i--) {
                if(Math.abs(corrs[n][logic.input.get(i)]) < GENE_CORREL_THRESH_MIN) {
                    logic.remove(i);
                }
            }
            if(logic.input.size() > 0) {
                logic.O_default.set(1 - logic.O.get(logic.input.size() - 1));
            }
            
            //update logic table of "n" again
            netD.nodes.get(n).update();
        }
    }
    
    private static String[] createFolders(String workDir, String rep) {
        String[] dirs = {DIR_OUTPUT, DIR_OUTPUT_RULE, DIR_OUTPUT_STRUCTS};

        for (int i = 0; i < dirs.length; i++) {
            dirs[i] = workDir + "\\" + dirs[i];
            dirs[i] = dirs[i].replace("out", rep);
            File theDir = new File(dirs[i]);

            // if the directory does not exist, create it
            if (!theDir.exists()) {
                try {
                    theDir.mkdir();
                    System.out.println("Created directory: " + theDir.getName());
                } catch (SecurityException se) {
                    //handle it
                }
            }
        }
        
        return dirs;
    }
    
    public ArrayList<NetInfo> GA_searchBN(int noStates, 
            int maxIteration, int maxConverge, int populationSize, double eliteRatio) throws InterruptedException {
        int numGenes = g_initial_data.numGenes;
//        boolean[][] data = g_initial_data.getData();
        ArrayList<NetInfo> population = new ArrayList<NetInfo>();
        int eliteSize = (int)Math.ceil(eliteRatio * populationSize);
        int noNet;
                
        for(noNet = 1; noNet <= populationSize; noNet++) {
            String netName = "RBN" + String.valueOf(noNet);
            NetInfo initNet = this.createRBN_BA(netName, numGenes, BA_EDGE_TO_ADD, 
                    BA_PROBABILITY, g_initial_data);                
            
            if(g_initial_data.finalSourceDest == false) {
                this.refineNetwork(initNet.netD);
            }
//            Utils.printInputs(initNet.netD.nodes);
            population.add(initNet);
        }
        Config.out("GA_searchBN", "Created an initial population of random boolean networks!");
        String stateSet = NetData.generateInitialStates("RBN1", String.valueOf(noStates));                        
                
                
        double bestFitness = Double.MIN_VALUE;        
        int noIter = 0;
        int noConvergeF = 0;
        boolean done = false;
        
        while (! done) {
            ++ noIter;                        
//            System.out.println("Processing the network: " + net.netD.networkName + " -------------------------------------------------------------------->>>");
//            System.out.println("Calculating scores .......");
            
            /*for(int i = 0; i < population.size(); i++) {
                NetInfo net = population.get(i);                
                if(net.fitness >= 0) continue;
                
                boolean succ = InferBN.calScore(net, stateSet, this.g_initial_data);
                
                if (succ == false || net.fitness < 0) {
                    Config.out("GA_searchBN", "Failed to calculate score of the network: " + net.netD.networkName);
                    population.remove(i);
                    -- i;                    
                }
            }*///non-parallel evaluation
            
            //parallel evaluation
            InferBN.paraEvaluate(population, stateSet, this.g_initial_data);
            
            for(int i = 0; i < population.size(); i++) {
                NetInfo net = population.get(i);                
                
                if (net.fitness < 0) {
                    Config.out("GA_searchBN", "Failed to calculate score of the network: " + net.netD.networkName);
                    population.remove(i);
                    -- i;                    
                }
            }
            //END: parallel evaluation
            
            if(population.isEmpty()) {
                Config.out("GA_searchBN", "Failed to calculate score of all networks! Stoped!");
                break;
            }
            
            Utils.quickSort(population, 0, population.size() - 1);
//            System.out.printf("Iter %d - scores = ", noIter);
//            for(int i = 0; i < population.size(); i++) {
//                System.out.print(population.get(i).fitness + " ");
//            }
//            System.out.println();            
            
            double cur_fitness = population.get(population.size() - 1).fitness;
            Config.out("GA_searchBN", "Iter " + noIter + " - best fitness = " + cur_fitness);
                        
            if(Math.abs(bestFitness - cur_fitness) <= SCORE_PRECISION) {
                ++ noConvergeF;
                
                if(noConvergeF >= maxConverge) {
                    Config.out("GA_searchBN", "The best fitness is converged! Stoped!");
                    break;
                }
            } else {
                if(cur_fitness - bestFitness > SCORE_PRECISION) {
                    bestFitness = cur_fitness;
                    noConvergeF = 1;                    
                } else {
                    //cur_fitness is bad                    
                }
            }
            
            if(Math.abs(cur_fitness - 1) <= SCORE_PRECISION) {
                Config.out("GA_searchBN", "Found perfect or nearly perfect networks! Stoped!");
                break;
            }                        
            
            if(noIter >= maxIteration) {
                Config.out("GA_searchBN", "The maximum number of iterations is reached! Stoped!");
                break;
            }                                    
            
            //Update number of convergences for each gene
            //not yet            
            
            //Create new population
            ArrayList<NetInfo> new_population = new ArrayList<NetInfo>();
            //Keep ELITE members of the current population
            if(population.size() <= eliteSize) {
                new_population.addAll(population);
            } else {
                for(int i = population.size() - eliteSize; i < population.size(); i++) {
                    new_population.add(population.get(i));
                }
            }
            
            //Mutate population
            int noMutations = populationSize - new_population.size();
            double sumFitness = Utils.sum(population);
            double[] probs_select = new double[population.size()];
            int[] index_select = new int[population.size()];            
            int cnt = 0;
            boolean[] nonExistNeighbor = new boolean[population.size()];
            
            for (int i = 0; i < population.size(); i++) {
                probs_select[i] = population.get(i).fitness / sumFitness;
                index_select[i] = i;
                nonExistNeighbor[i] = false;
            }
            EnumeratedIntegerDistribution rand = new EnumeratedIntegerDistribution(rand_select_member, 
                    index_select, probs_select);
        
//            System.out.println("Creating new population .......");
            while(cnt < noMutations) {
                int iNet = rand.sample();
                if(nonExistNeighbor[iNet] == true) continue;
                
                ArrayList<NetInfo> neighborNets = this.getNeighborNetworks(population.get(iNet), 
                        1, noNet + 1);
                
                if(neighborNets.size() > 0) {
                    for (int i = 0; i < neighborNets.size(); i++) {
                        new_population.add(neighborNets.get(i));
                        Calc.datas.put(neighborNets.get(i).netD.networkName, neighborNets.get(i).netD);
                    }
                    noNet += neighborNets.size();
                    cnt += neighborNets.size();
                    
                } else {
                    nonExistNeighbor[iNet] = true;
                }
                
                if(! Utils.exist(false, nonExistNeighbor)) break;
            }
            if(cnt == 0) {
                Config.out("GA_searchBN", "Cannot produce new members from the current population! Stoped!");
                break;
            }
                                                                                    
            //clear the current population
            //Remove the networks from Calc, except ELITE networks
            for(int i = 0; i < population.size() - eliteSize; i++) {
                Calc.datas.remove(population.get(i).netD.networkName);
                population.get(i).free();                
            }
            population.clear();
            
            population = new_population;
            System.gc();
        }
        
        //select optimal networks
        ArrayList<NetInfo> optimalNets = new ArrayList<NetInfo>();
        
        if(! population.isEmpty()) {
            bestFitness = population.get(population.size() - 1).fitness;
            optimalNets.add(population.get(population.size() - 1));
            
            for (int i = population.size() - 2; i >= 0; i--) {
                if(Math.abs(bestFitness - population.get(i).fitness) <= SCORE_PRECISION) {
                    optimalNets.add(population.get(i));
                } else {
                    break;
                }
            }
            
            //clear the current population            
            for (int i = 0; i < population.size() - optimalNets.size(); i++) {
                Calc.datas.remove(population.get(i).netD.networkName);
                population.get(i).free();
            }
            population.clear();
        }
        population = null;
        System.gc();                
        
        Config.out("GA_searchBN", "Found " + optimalNets.size() + " optimal networks with fitness = " + bestFitness);
        return optimalNets;
    }
        
    public static void paraEvaluate(ArrayList<NetInfo> population, 
            String stateSet, GeneBData g_initial_data) throws InterruptedException {
        
        ArrayList<ComputeUnit> runnables = new ArrayList<ComputeUnit>();
        
        for (int i = 0; i < population.size(); i++) {
            NetInfo net = population.get(i);
            if (net.fitness >= 0) continue;
            
            runnables.add(new ComputeUnit(net, stateSet, g_initial_data));
        }

        int threads = runnables.size();
        ExecutorService exec = Executors.newFixedThreadPool(threads);
        CountDownLatch latch = new CountDownLatch(threads);
        
        for (ComputeUnit r : runnables) {
            r.setLatch(latch);
            exec.execute(r);
        }

        latch.await();
        exec.shutdown();        
    }
    
    //Random generators    
    Random rand_select_tarGene;
    RandomGenerator rand_select_input;
    RandomGenerator rand_select_member;
    
    EnumeratedIntegerDistribution rand_pert_type;    
    int[] perts_ToGenerate = new int[]{0, 1, 2, 3, 4, 5};
    double[] perts_Prob = new double[]{0.18, 0.18, 0.18, 0.18, 0.18, 0.1};
    
    GeneBData g_initial_data;
    int[][] fixedInputs;
        //each row contains a list of fixed regulators for each target gene        
    boolean[] fixedGenes;
    int gScore_thresh;
    double[] prev_g_scores;
    int[] g_noConverges;
    
    Set<String> ExaminingRules;
}
