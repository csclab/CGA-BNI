/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package bni;

import bni.comp.GeneBData;
import bni.comp.NetInfo;
import bni.comp.Regulator;
import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.Iterator;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.TreeSet;
import mod.jmut.core.comp.Interaction;
import mod.jmut.core.comp.LogicTable;
import mod.jmut.core.comp.NetData;
import mod.jmut.core.comp.Node;
import org.apache.commons.math3.distribution.TDistribution;
import org.apache.commons.math3.stat.correlation.PearsonsCorrelation;

/**
 *
 * @author colin
 */
public class Utils {
    public static String[][] loadTextFile(String filePath, String delim, boolean ignoreFirstLine) {
        String [][] data;   //foramt: [col][row]
        int iL = 0;
        try {
            int numLines = countLines(filePath);
            int numCols = countColumns(filePath, delim);
//            System.out.println("numCols=" + numCols);
            
            if(ignoreFirstLine) {
                numLines = numLines - 1;
            }
            data = new String[numCols][numLines];
            
            BufferedReader input = new BufferedReader(new FileReader(filePath));
            String line;

            if(ignoreFirstLine) {
                input.readLine();   //ignore first line (title line)
            }
                        
            while ((line = input.readLine()) != null) {
                if (line.length() == 0) {
                    continue;
                }

                String[] params = line.split(delim, -1);
                for(int c = 0; c < params.length; c ++) {
                    data[c][iL] = params[c];
                }
                ++ iL;
            }
            
            input.close();
            
        } catch (Exception ex) {
            ex.printStackTrace();
            System.out.println("iL=" + iL);
            return null;
        }
        return data;
    }
    
    public static String[][] loadTextFile_byRow(String filePath, String delim, boolean ignoreFirstLine) {
        String [][] data;   //foramt: [row][col]
        int iL = 0;
        try {
            int numLines = countLines(filePath);            
            
            if(ignoreFirstLine) {
                numLines = numLines - 1;
            }
            data = new String[numLines][];
            
            BufferedReader input = new BufferedReader(new FileReader(filePath));
            String line;

            if(ignoreFirstLine) {
                input.readLine();   //ignore first line (title line)
            }
                        
            while ((line = input.readLine()) != null) {
                if (line.length() == 0) {
                    continue;
                }

                String[] params = line.split(delim, -1);
                data[iL] = new String[params.length];
                
                for(int c = 0; c < params.length; c ++) {
                    data[iL][c] = params[c];
                }
                ++ iL;
            }
            
            input.close();
            
        } catch (Exception ex) {
            ex.printStackTrace();
            System.out.println("iL=" + iL);
            return null;
        }
        return data;
    }
    
    private static int countLines(String filePath) {
        int count = 0;
        
        try {
            BufferedReader input = new BufferedReader(new FileReader(filePath));
            String line;            
            
            while ((line = input.readLine()) != null) {
                if (line.length() > 0) {
                    ++ count;
                }                                
            }
            
            input.close();
            
        } catch (Exception ex) {
            ex.printStackTrace();
            return 0;
        }
        return count;
    }
    
    private static int countColumns(String filePath, String delim) {
        int count = 0;
        
        try {
            BufferedReader input = new BufferedReader(new FileReader(filePath));
            String line;            
            
            input.readLine();
            while ((line = input.readLine()) != null) {
                if (line.length() > 0) {
                    String[] params = line.split(delim);                    
                    count = params.length;
                    break;
                }                                
            }
            
            input.close();
            
        } catch (Exception ex) {
            ex.printStackTrace();
            return 0;
        } 
        return count;
    }
    
    public static GeneBData loadOriginNetwork(String path, int numGenes) {
        String [][] data = Utils.loadTextFile(path, "\\t", false);
        GeneBData gdata = new GeneBData(numGenes);
        
        gdata.loadGoldData(data);
        return gdata;
    }
    
    //Stastic
    public static boolean[][] parseBoolData(String [][] data, int colStart, int colEnd) {
        int numLines = data[0].length;
        int numCols = colEnd - colStart + 1;
        
        boolean[][] results = new boolean[numLines][numCols];
        
        for(int l = 0; l < numLines; l ++) {
            for(int c = colStart; c <= colEnd; c ++) {
                results[l][c - colStart] = data[c][l].equals("1")? true: false;
            }
        }
        return results;
    }
    
    public static int[] parseIntData(String [][] data, int col) {
        int numLines = data[0].length;        
        
        int[] results = new int[numLines];
        
        for(int l = 0; l < numLines; l ++) {            
            results[l] = Integer.parseInt(data[col][l]);
        }
        return results;
    }
    
    public static double[] parseDoubleData(String [][] data, int col) {
        int numLines = data[0].length;        
        
        double[] results = new double[numLines];
        
        for(int l = 0; l < numLines; l ++) {            
            results[l] = Double.parseDouble(data[col][l]);
        }
        return results;
    }
    
    public static boolean[] parse(boolean[][] data, int colIndex) {
        int numLines = data.length;
        
        boolean [] results = new boolean[numLines];
        for(int l = 0; l < numLines; l ++) {
            results[l] = data[l][colIndex];
        }
        return results;
    }
    
    public static double[][] toDouble(boolean[][] data) {
        int numLines = data.length;
        int numCols = data[0].length;
        
        double[][] results = new double[numLines][numCols];
        
        for(int l = 0; l < numLines; l ++) {
            for(int c = 0; c < numCols; c ++) {
                if(data[l][c] == true) results[l][c] = 1;
                else results[l][c] = 0;
            }
        }
        return results;
    }
    
    public static int[] makeIntArray(int len, int dvalue) {
        int[] arr = new int[len];
        
        for (int i = 0; i < arr.length; i++) {
            arr[i] = dvalue;
        }
        return arr;
    }
    
    public static double[] calCorrelation(double[] v1, double[] v2, int size, int numTail) {
        PearsonsCorrelation corr = new PearsonsCorrelation();        
        double coeff = corr.correlation(v1, v2);
        double t = coeff * Math.sqrt((size - 2) / (1 - coeff * coeff));
        
	double p = 1 - new TDistribution(size - 2).cumulativeProbability(Math.abs(t));
        p = numTail * p;
        
        return new double[]{coeff, p};
    }
    
    public static int minus(boolean b1, boolean b2) {
        if(b1 == b2) return 0;
        else {
            if(b1 == true) return 1;
            else return -1;
        }
    }
    
    public static double sum(double[] someArray) {
        double sum = 0;

        for (double i : someArray) {
            sum += i;
        }

        return sum;
    }
    
    public static double mean(double[] someArray) {
        return Utils.sum(someArray) / someArray.length;
    }
    
    public static double std(double[] someArray, double mean) {
        double res = 0;
        int N = someArray.length;
        
        for (int i = 0; i < N; i++) {
            res += Math.pow(someArray[i] - mean, 2);
        }
        
        res = Math.sqrt(res / N);
        return res;
    }
    
    public static double[] toArray(List<Double> doubles) {
        double[] target = new double[doubles.size()];
        
        for (int i = 0; i < target.length; i++) {            
            target[i] = doubles.get(i);
        }
        return target;
    }
    
    public static double sum(ArrayList<NetInfo> someArray) {
        double sum = 0;

        for (NetInfo net : someArray) {
            sum += net.fitness;
        }

        return sum;
    }
    
    public static boolean exist(int value, int[] array) {
        boolean found = false;
        
        for(int v: array) {
            if(value == v) {
                found = true;
                break;
            }
        }
        return found;
    }
    
    public static boolean exist(boolean value, boolean[] array) {
        boolean found = false;
        
        for(boolean v: array) {
            if(value == v) {
                found = true;
                break;
            }
        }
        return found;
    }
    
    public static boolean exist(int value, ArrayList<Integer> array) {
        boolean found = false;
        
        for(int v: array) {
            if(value == v) {
                found = true;
                break;
            }
        }
        return found;
    }
    
    public static boolean exist(int gene, Set<Regulator> regulators) {
        boolean found = false;
        
        for(Regulator rg: regulators) {
            if(gene == rg.regulator) {
                found = true;
                break;
            }
        }
        return found;
    }
    
    public static int index(int value, int[] array) {
        int id = -1;
        
        for (int i = 0; i < array.length; i++) {
            if(value == array[i]) {
                id = i;
                break;
            }
        }
        return id;
    }
    
    public static int index(String value, String[] array) {
        int id = -1;
        
        for (int i = 0; i < array.length; i++) {
            if(value.equalsIgnoreCase(array[i])) {
                id = i;
                break;
            }
        }
        return id;
    }
    
    public static int[] convertIntegers(List<Integer> integers) {
        int[] ret = new int[integers.size()];
        for (int i = 0; i < ret.length; i++) {
            ret[i] = integers.get(i).intValue();
        }
        return ret;
    }
    
    public static boolean[][] convert_2DBool(String state) {
        int numNodes = state.length();
        boolean[][] bstate = new boolean[1][numNodes];
        
        for(int i = 0; i < numNodes; i ++) {
            bstate[0][i] = state.charAt(i) == '0' ? false : true;
        }    
        return bstate;
    }
    
    public static int countFixedRegulators(TreeSet<Regulator> regGenes) {
        int cnt = 0;

        for (Regulator reg : regGenes) {
            if(reg.rank == Regulator.RANK_DIRECT_POSITIVE
                        || reg.rank == Regulator.RANK_DIRECT_NEGATIVE) {
                ++cnt;
            }
        }

        return cnt;
    }
    
    public static void printFixedInputs(int[][] fixedInputs) {
        int numGenes = fixedInputs.length;
        
        System.out.println("Fixed inputs of each gene:");
        
        for (int g = 0; g < numGenes; g++) {
            System.out.printf("Gene %d:\t", g + 1);
            
            for(int i = 0; i < fixedInputs[g].length; i++) {
                System.out.printf("%d ", fixedInputs[g][i] + 1);
            }
            System.out.println();
        }
    }
    
    public static void printInputs(ArrayList<Node> nodes) {
        for(int n = 0; n < nodes.size(); n++) {
            LogicTable logic = nodes.get(n).getLogicTable();
            int noInputs = logic.input.size();
            
            System.out.printf("Gene %d:\t", n + 1);
            
            for(int i = 0; i < noInputs; i++) {
                System.out.printf("%d ", logic.input.get(i) + 1);
            }
            System.out.println();
        }
    }
    
    public static void outputDynamicsAccuracy(String path, ArrayList<NetInfo> optimalNets, long searchTime, String del) throws FileNotFoundException {
        PrintWriter pw_net = new PrintWriter(new FileOutputStream(path + "_net.csv"),true);                
        PrintWriter pw_node = new PrintWriter(new FileOutputStream(path + "_node.csv"),true);                
        
        pw_net.println("Search time" + del + searchTime);
        pw_net.println("Network" + del + "Dynamics Accuracy");        
        pw_node.println("Network" + del + "Gene" + del + "Absolute distance");        
        
        for(int i = 0; i < optimalNets.size(); i ++) {
            NetInfo net = optimalNets.get(i);
            pw_net.println(net.netD.networkName + del + (1 - net.avg_score));
            
            for(int g = 0; g < net.netD.nodes.size(); g ++) {
                pw_node.println(net.netD.networkName + del + net.netD.nodes.get(g).NodeID
                        + del + net.g_scores[g]);
            }
        }
        
        pw_net.close();
        pw_node.close();
    }
    
    public static void outputRules(String dir, NetData netD) {
        
        String filename = dir + netD.networkName + "_rules_0";
        Node.outputRules(filename, netD.nodes);
    }
    
    public static void updateNodes(NetData newNetD) {
        for (int n = 0; n < newNetD.nodes.size(); n++) {
            newNetD.nodes.get(n).update();
        }
    }
    
    public static void copy(int[][] sa, int[][] da) {
        for (int i = 0; i < sa.length; i++) {
            System.arraycopy(sa[i], 0, da[i], 0, sa[i].length);
        }
    }
    
    public static boolean[][] copy(boolean[][] origin) {
        boolean[][] cop = new boolean[origin.length][origin[0].length];

        for (int i = 0; i < cop.length; i++) {
            System.arraycopy(origin[i], 0, cop[i], 0, cop[i].length);
        }       
        
        return cop;
    }
    
    public static int[] copy(int[] origin) {
        int[] cop = new int[origin.length];
        
        System.arraycopy(origin, 0, cop, 0, cop.length);
        
        return cop;
    }
    
    public static boolean modifyLink(int a, int type, int b, int[][] C, 
                               int[][] posPaths, int[][] negPaths, boolean force, int mode) {
        int numGenes = C.length;
        int[][] new_posPaths = new int[numGenes][numGenes];
        int[][] new_negPaths = new int[numGenes][numGenes];
        
        Utils.copy(posPaths, new_posPaths);
        Utils.copy(negPaths, new_negPaths);
        
        //modify
        for (int s = 0; s < numGenes; s++) {
            for (int d = 0; d < numGenes; d++) {
                if(s == d || s == b || d == a) continue;
                
                int[][] arr1 = new_posPaths;
                int[][] arr2 = new_negPaths;
                if(type == -1) {
                    arr1 = new_negPaths;
                    arr2 = new_posPaths;
                }
                
                if (s == a && d == b) {
                    switch(mode) {
                        case Regulator.MODE_ADD:
                            arr1[s][d] += 1;
                            break;
                            
                        case Regulator.MODE_DEL:
                            arr1[s][d] -= 1;
                            break;
                            
                        case Regulator.MODE_CHG://pos to neg
                            arr1[s][d] -= 1;
                            arr2[s][d] += 1;
                            break;
                            
                        default:
                            break;
                    }                    
                } else {
                    if (s == a) {
                        switch (mode) {
                            case Regulator.MODE_ADD:
                                arr1[s][d] += posPaths[b][d];
                                arr2[s][d] += negPaths[b][d];
                                break;

                            case Regulator.MODE_DEL:
                                arr1[s][d] -= posPaths[b][d];
                                arr2[s][d] -= negPaths[b][d];
                                break;

                            case Regulator.MODE_CHG:
                                arr1[s][d] -= posPaths[b][d];
                                arr2[s][d] -= negPaths[b][d];
                                
                                arr1[s][d] += negPaths[b][d];
                                arr2[s][d] += posPaths[b][d];
                                break;

                            default:
                                break;
                        }                        
                    } else {
                        if (d == b) {
                            switch (mode) {
                                case Regulator.MODE_ADD:
                                    arr1[s][d] += posPaths[s][a];
                                    arr2[s][d] += negPaths[s][a];
                                    break;

                                case Regulator.MODE_DEL:
                                    arr1[s][d] -= posPaths[s][a];
                                    arr2[s][d] -= negPaths[s][a];
                                    break;

                                case Regulator.MODE_CHG:
                                    arr1[s][d] -= posPaths[s][a];
                                    arr2[s][d] -= negPaths[s][a];

                                    arr1[s][d] += negPaths[s][a];
                                    arr2[s][d] += posPaths[s][a];
                                    break;

                                default:
                                    break;
                            }
                        } else {
                            switch (mode) {
                                case Regulator.MODE_ADD:
                                    arr1[s][d] += posPaths[s][a] * posPaths[b][d] + negPaths[s][a] * negPaths[b][d];
                                    arr2[s][d] += posPaths[s][a] * negPaths[b][d] + negPaths[s][a] * posPaths[b][d];
                                    break;

                                case Regulator.MODE_DEL:
                                    arr1[s][d] -= posPaths[s][a] * posPaths[b][d] + negPaths[s][a] * negPaths[b][d];
                                    arr2[s][d] -= posPaths[s][a] * negPaths[b][d] + negPaths[s][a] * posPaths[b][d];
                                    break;

                                case Regulator.MODE_CHG://pos to neg
                                    arr1[s][d] -= posPaths[s][a] * posPaths[b][d] + negPaths[s][a] * negPaths[b][d];
                                    arr2[s][d] -= posPaths[s][a] * negPaths[b][d] + negPaths[s][a] * posPaths[b][d];
                                    
                                    arr1[s][d] += posPaths[s][a] * negPaths[b][d] + negPaths[s][a] * posPaths[b][d];
                                    arr2[s][d] += posPaths[s][a] * posPaths[b][d] + negPaths[s][a] * negPaths[b][d];
                                    break;

                                default:
                                    break;
                            }
                        }
                    }
                }                                                     
            }
        }
        
        if(force == true) {
            Utils.copy(new_posPaths, posPaths);
            Utils.copy(new_negPaths, negPaths);
            return true;
        }
        
        for (int s = 0; s < numGenes; s++) {
            for (int d = 0; d < numGenes; d++) {
                if(s == d || s == b || d == a) continue;
                
                if(new_posPaths[s][d] != posPaths[s][d]
                   || new_negPaths[s][d] != negPaths[s][d]) {
                    int old_dist = Utils.pairOK(s, d, C, posPaths, negPaths);
                    int new_dist = Utils.pairOK(s, d, C, new_posPaths, new_negPaths);
                    
                    if(new_dist == -1) {    //OKE
                        continue;
                    } else {
                        if(old_dist != -1 && new_dist <= old_dist) {
                            //Both old & new results are bad, but the new one has an improvement
                            continue;
                        } else {
                            return false;
                        }
                    }
                }
            }
        }
        
        Utils.copy(new_posPaths, posPaths);
        Utils.copy(new_negPaths, negPaths);
        return true;
    }
    
    public static int pairOK(int s, int d, int[][] C, int[][] posPaths, int[][] negPaths) {
        int dist = -1;
        
        switch (C[s][d]) {
            case Regulator.RANK_INDIRECT_POSITIVE:
                dist = posPaths[s][d] - negPaths[s][d];
                if(dist > 0) {
                    dist = -1;  //OKE
                } else {
                    dist = Math.abs(dist);
                }
                break;
                
            case Regulator.RANK_INDIRECT_NEGATIVE:
                dist = posPaths[s][d] - negPaths[s][d];
                if(dist < 0) {
                    dist = -1;  //OKE
                }
                break;
                
            case Regulator.RANK_INDIRECT_NEUTRAL:
                dist = Math.abs(posPaths[s][d] - negPaths[s][d]);
                if(dist == 0) {
                    dist = -1;  //OKE
                }
                break;
            
            case Regulator.RANK_NOLINK:
                if(posPaths[s][d] == 0 && negPaths[s][d] == 0) {
                    dist = -1;  //OKE
                } else {
                    dist = posPaths[s][d] + negPaths[s][d];
                }
                break;
                
            default:
                break;
        }
        
        return dist;
    }
    
    public static void quickSort(ArrayList<NetInfo> A, int lower, int upper){                
        double x = A.get((lower + upper) / 2).fitness;
        int i = lower;
        int j = upper;
        
        while(i <= j){
            while(Double.compare(A.get(i).fitness, x) < 0) i++;
            while(Double.compare(A.get(j).fitness, x) > 0) j--;
            
            if (i <= j){                
                NetInfo temp = A.get(i);
                A.set(i, A.get(j));
                A.set(j, temp);

                i++;
                j--;
            }
            //System.out.println("i = " + i + ", j = " + j);
        }
        
        if (j > lower) quickSort(A, lower, j);
        if (i < upper) quickSort(A, i, upper);
    }    
    
    /*
     * Other tools
     */
    public static GeneBData load_GENEI3_results(String path, int numGenes, int maxNoEdges) {
        String [][] data = Utils.loadTextFile(path, "\\t", false);
        GeneBData gdata = new GeneBData(numGenes);
        
        int numLines = data[0].length;
        double[] scores = Utils.parseDoubleData(data, 2);
        double mean = Utils.sum(scores) / numLines;
        System.out.println("[load_GENEI3_results] mean = " + mean);
        
        HashMap<Integer, TreeSet<Regulator>> regulators = gdata.getRegulators();
        if(maxNoEdges > numLines) maxNoEdges = numLines;
        
        for(int l = 0; l < maxNoEdges; l ++) {
            if(scores[l] < mean) {
                break;
            }
            
            int srcGene = Integer.parseInt(data[0][l].substring(1));
            int tarGene = Integer.parseInt(data[1][l].substring(1));
            int type = 1;            
                        
            regulators.get(tarGene).add(new Regulator(srcGene, type, 0));
        }
        
        gdata.finalSourceDest = true;
        return gdata;                        
    }
    
    public static GeneBData load_ARACNE_results(String path, int numGenes, int maxNoEdges) {
        String [][] data = Utils.loadTextFile_byRow(path, "\\t", false);
        GeneBData gdata = new GeneBData(numGenes);
        
        int numLines = data.length;
        int type = 1; 
//        double[] scores = Utils.parseDoubleData(data, 2);
//        double mean = Utils.sum(scores) / numLines;
//        System.out.println("[load_GENEI3_results] mean = " + mean);
                
        HashMap<String, Double> passedMap = new HashMap<String, Double>();
        
        for(int li = 0; li < numLines; li ++) {

            if(data[li].length < 3) continue;
            
            int srcGene = Integer.parseInt(data[li][0].substring(1));
            
            for(int tg = 1; tg < data[li].length; tg += 2) {
                int tarGene = Integer.parseInt(data[li][tg].substring(1));                           

                double mi = Double.parseDouble(data[li][tg + 1]);
                passedMap.put(srcGene + "_" + tarGene, mi);
            }
        }
        
        LinkedHashMap<String, Double> sortedMap = sortHashMapByValues(passedMap);
        HashMap<Integer, TreeSet<Regulator>> regulators = gdata.getRegulators();   
        int cnt = 0;
        
        for (Map.Entry<String, Double> entry : sortedMap.entrySet()) {
            String key = entry.getKey();
            System.out.println("[ARACNE] " + key + "\t" + entry.getValue());
            
            String[] params = key.split("_|\\.", -1);
            int srcGene = Integer.parseInt(params[0]);
            int tarGene = Integer.parseInt(params[1]);
            
            regulators.get(tarGene).add(new Regulator(srcGene, type, 0));
            
            ++ cnt;
            if(cnt >= maxNoEdges) break;
        }                             
        
        gdata.finalSourceDest = true;
        return gdata;                        
    }
    
    public static GeneBData load_BC3NET_results(String path, int numGenes) {
        String [][] data = Utils.loadTextFile_byRow(path, "\\s+", false);
        GeneBData gdata = new GeneBData(numGenes);
        
        int numLines = data.length;
        int type = 1; 
//        double[] scores = Utils.parseDoubleData(data, 2);
//        double mean = Utils.sum(scores) / numLines;
//        System.out.println("[load_GENEI3_results] mean = " + mean);
        
        HashMap<Integer, TreeSet<Regulator>> regulators = gdata.getRegulators();        
        
        for(int li = 0; li < numLines; li ++) {
            if(data[li].length < 3) continue;
            
            int srcGene = Integer.parseInt(data[li][0].substring(1));                        
            int tarGene = Integer.parseInt(data[li][1].substring(1));                           

            regulators.get(tarGene).add(new Regulator(srcGene, type, 0));
            regulators.get(srcGene).add(new Regulator(tarGene, type, 0));
        }
        
        gdata.finalSourceDest = true;
        return gdata;                        
    }
    
    public static boolean checkExistInteraction(Interaction newina, int numofina, ArrayList<Interaction> ina) {
        boolean exist = false;
        int i;
        for (i = 0; i < numofina; i++) {
            if (newina.NodeSrc.compareTo(ina.get(i).NodeSrc) == 0 && newina.NodeDst.compareTo(ina.get(i).NodeDst) == 0) {
                exist = true;
                break;
            }
        }
        return exist;
    }

    private static String toCSVString(String s, String del) {
        int len = s.length();
        StringBuilder sb = new StringBuilder();
        
        for(int i = 0; i < len - 1; i++) {
            sb.append(s.charAt(i)).append(del);
        }
        sb.append(s.charAt(len - 1));
        
        return sb.toString();
    }
    
    public static void outputKOData(ArrayList<String> wtAtts, String[] pertAtts, int[] wtIndexs,
            String networkName, String path, String del) throws FileNotFoundException {
        
        String filePath = path + networkName + "_knockouts.tsv.bool.csv";
        PrintWriter pw_net = new PrintWriter(new FileOutputStream(filePath),true);                
        
        for(int i = 0; i < wtAtts.size(); i++) {
            pw_net.println("-1" + del + i + del + toCSVString(wtAtts.get(i), del));
        }
        
        for(int i = 0; i < pertAtts.length; i++) {
            pw_net.println(String.valueOf(i + 1) + del + wtIndexs[i] + del + 
                                toCSVString(pertAtts[i], del));
        }
        
        pw_net.close();
    }

    public static void outputGoldNetwork(NetData netD,
            String path, String del) throws FileNotFoundException {
        
        String filePath = path + netD.networkName + "_goldstandard_signed.tsv";
        PrintWriter pw_net = new PrintWriter(new FileOutputStream(filePath),true);                
        
        for(int i = 0; i < netD.rndina.size(); i++) {
            Interaction ina = netD.rndina.get(i);
            String src = "G" + String.valueOf( Integer.valueOf(ina.NodeSrc) + 1 );
            String dst = "G" + String.valueOf( Integer.valueOf(ina.NodeDst) + 1 );
            String sign = "+";
            if(ina.InteractionType == -1) sign = "-";
            
            pw_net.println(src + del + dst + del + sign);
        }
        pw_net.close();
    }
    
    public static LinkedHashMap<String, Double> sortHashMapByValues(
            HashMap<String, Double> passedMap) {
        List<String> mapKeys = new ArrayList(passedMap.keySet());
        List<Double> mapValues = new ArrayList(passedMap.values());
        
        Collections.sort(mapValues, Collections.reverseOrder());
        Collections.sort(mapKeys);

        LinkedHashMap<String, Double> sortedMap =
                new LinkedHashMap();

        Iterator<Double> valueIt = mapValues.iterator();
        while (valueIt.hasNext()) {
            double val = valueIt.next();
            Iterator<String> keyIt = mapKeys.iterator();

            while (keyIt.hasNext()) {
                String key = keyIt.next();
                Double comp1 = passedMap.get(key);
                Double comp2 = val;

                if (comp1.equals(comp2)) {
                    keyIt.remove();
                    sortedMap.put(key, val);
                    break;
                }
            }
        }
        return sortedMap;
    }
}
