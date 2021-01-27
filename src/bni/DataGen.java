/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package bni;

import bni.comp.GeneBData;
import bni.comp.NetInfo;
import java.awt.Point;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Random;
import mod.Knockout;
import mod.Mutation;
import mod.OverExpression;
import mod.jmut.core.Calc;
import mod.jmut.core.Config;
import mod.jmut.core.Util;
import mod.jmut.core.comp.Attractor;
import mod.jmut.core.comp.Interaction;
import mod.jmut.core.comp.NetData;

/**
 *
 * @author colin_PC2
 */
public class DataGen {
    
    public static NetData createRBN_BA(int numGenes, int edgesToAdd, int numofinitnodes, 
            double probability, String netName) {
        ArrayList<Interaction> inatemp = new ArrayList<Interaction>();
        
        try {
            Random random = new Random();
            int i, j;
            int degrees[] = new int[numGenes];

            int numofedges = 0;                        
            int inatype;            

            for (i = 0; i < numofinitnodes; i++) {
                for (j = (i + 1); j < numofinitnodes; j++) {
                    inatype = (Math.random() < probability) ? -1 : 1;
                    inatemp.add(numofedges, new Interaction());

                    inatemp.get(numofedges).InteractionType = inatype;
                    if (Math.random() < 0.5) {
                        inatemp.get(numofedges).NodeSrc = Integer.toString(i);
                        inatemp.get(numofedges).NodeDst = Integer.toString(j);
                    } else {
                        inatemp.get(numofedges).NodeDst = Integer.toString(i);
                        inatemp.get(numofedges).NodeSrc = Integer.toString(j);

                    }
                    degrees[i]++;
                    degrees[j]++;
                    numofedges++;
                }
            }
            
            for (i = numofinitnodes; i < numGenes; i++) {
                double degreeIgnore = 0;
                double oldTotalDegrees = 2.0d * numofedges;
                
                for (int m = 0; m < edgesToAdd; m++) {
                    double prob = 0;
                    double randNum = random.nextDouble();
                    for (j = 0; j < i; j++) {
                        boolean existing = true;
                        Interaction temp;
                        int numoftrial = 0;
                        while (true) {
                            numoftrial++;
                            inatype = (Math.random() < probability) ? -1 : 1;

                            temp = new Interaction();

                            temp.InteractionType = inatype;
                            if (Math.random() < 0.5) {
                                temp.NodeSrc = Integer.toString(i);
                                temp.NodeDst = Integer.toString(j);
                            } else {
                                temp.NodeDst = Integer.toString(i);
                                temp.NodeSrc = Integer.toString(j);

                            }

                            if (Utils.checkExistInteraction(temp, numofedges, inatemp) == false) {
                                prob += (double) ((double) degrees[j]) / ((double) (oldTotalDegrees/*2.0d * numofedges*/) - degreeIgnore);
                                existing = false;
                                break;
                            }
                            if (edgesToAdd >= 2 && numoftrial > 10) {//colin: fix hang out bug
                                break;
                            }
                        }
                        if (randNum <= prob && existing == false) {
                            inatemp.add(numofedges, temp);

                            degreeIgnore += degrees[j];
                            degrees[i]++;
                            degrees[j]++;

                            numofedges++;

                            break;
                        }
                    }
                }
            }

        } catch (Exception e) {
            e.printStackTrace();
        }

        NetData.loadNetwork_v2(netName, inatemp, numGenes);
        NetData.generateCurrentRule(netName, 1);// OR rule only
        
        return Calc.datas.get(netName);        
    }

    public static boolean generate_KO_Attractors(NetData netD, int noStates, String path) throws FileNotFoundException {
        
        Calc cal = new Calc();
        String stateSet = NetData.generateInitialStates(netD.networkName, String.valueOf(noStates));        
        int numGenes = netD.nodes.size();        
        boolean[][] states = Calc.states.get(stateSet).getCurrent();
        
        //Find all wild-type attractors
        ArrayList<Attractor> AllAttractors = cal.findAttractors_v2(netD, states, false);
        
        if(AllAttractors == null) {
            Config.out("generate_KO_Attractors", "Failed to find all wild-type attractors of the network!");
            return false;
        }
        System.out.println("Number of wild-type atts = " + AllAttractors.size());
        
        ArrayList<String> wtAtts = new ArrayList<String>();
        String[] pertAtts = new String[numGenes];
        int[] wtIndexs = new int[numGenes];        
        boolean[] dones = new boolean[numGenes];
        
        for (int pertG = 0; pertG < numGenes; pertG++) {
            dones[pertG] = false;
        }
        
        for (int i = 0; i < AllAttractors.size(); i++) {
            if(! Utils.exist(false, dones)) break;
        
            Attractor att = AllAttractors.get(i);
            if (att.Length > 1) {
                continue;
            }
            
            boolean[][] wt_state = Utils.convert_2DBool(att.States.get(0));            
            boolean exist_KOFixedPointAtt = false;
            
            for (int pertG = 0; pertG < numGenes; pertG++) {
                if(dones[pertG] == true) continue;
                
                Mutation mut = new Knockout();

                Util.perturbNode(netD.nodes.get(pertG), mut);
                ArrayList<Attractor> newAttractors = cal.findAttractors_v2(netD, wt_state, false);
                if (newAttractors == null) {
                    Config.out("generate_KO_Attractors", "Failed to find the perturbed attractor along with the mutation on gene " + pertG);
                    continue;
                }

                netD.nodes.get(pertG).restoreOriginalLogicTable();

                //Add perturbed attractor
                Attractor patt = newAttractors.get(0);
                if(patt.Length == 1) {
                    pertAtts[pertG] = patt.States.get(0);
                    wtIndexs[pertG] = wtAtts.size();
                    dones[pertG] = true;
                    
                    exist_KOFixedPointAtt = true;
                }
            }
            
            if(exist_KOFixedPointAtt) {
                wtAtts.add(att.States.get(0));
            }
        }
        
        Calc.states.remove(stateSet);
        Calc.datas.remove(netD.networkName);
        
        if(! Utils.exist(false, dones)) {
            Utils.outputGoldNetwork(netD, path, "\t");
            Utils.outputKOData(wtAtts, pertAtts, wtIndexs, netD.networkName, path, ",");
            
            Config.out("generate_KO_Attractors", "Successfully output KO data and Gold network: " + netD.networkName);
            return true;
        } else {
            return false;
        }
    }
    
    public static final int START_INDEX_NETWORK = 30;
    
    public static void createRBNs(int numNodes, int numNetworks) throws FileNotFoundException, IOException, InterruptedException {
        
        File theDir = new File("Size" + numNodes + "\\");

        if (!theDir.exists()) {
            try {
                theDir.mkdir();
                System.out.println("Created directory: " + theDir.getName());
            } catch (SecurityException se) {
            }
        }

        int cnt = 0;
        while( cnt < numNetworks ) {
            String netName = "Size" + numNodes + "_RBN" + String.valueOf(cnt + START_INDEX_NETWORK);
            NetData netD = createRBN_BA(numNodes, 2, 3, 0.5, netName);
            
            boolean succ = generate_KO_Attractors(netD, InferBN.SCORE_NO_INITIAL_STATES, 
                    theDir.getName() + "\\");
            if(succ) {
                ++ cnt;
            }
        }
    }  
    
    public static void main(String[] args) throws FileNotFoundException, IOException, InterruptedException {
        int NO_SIZES = 10;
        int NO_NETWORKS = 20;
        
        int[] sizes = new int[NO_SIZES];
        
        for(int i = 0; i < NO_SIZES; i++) {
            sizes[i] = 10 * (i + 1);
            
            createRBNs(sizes[i], NO_NETWORKS);
        }
    }    

}
