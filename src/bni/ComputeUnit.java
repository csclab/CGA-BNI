/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package bni;

import bni.comp.GeneBData;
import bni.comp.NetInfo;
import java.util.concurrent.CountDownLatch;
import mod.jmut.core.Calc;

/**
 *
 * @author colin_PC2
 */
public class ComputeUnit implements Runnable {

    NetInfo net;
    boolean[][] states;
    boolean[][] base_data;
    int[] expGenes;
    int[] wildTypes;
            
    private CountDownLatch latch;

    public ComputeUnit(NetInfo net, String stateSet, GeneBData g_initial_data) {
        this.net = net;
        
        boolean[][] origin_states = Calc.states.get(stateSet).getCurrent();
        this.states = Utils.copy(origin_states);
        
        boolean[][] __base_data = g_initial_data.getData();
        int[] __expGenes = g_initial_data.getExpGenes();
        int[] __wildTypes = g_initial_data.getWildTypes();
        this.base_data = Utils.copy(__base_data);
        this.expGenes = Utils.copy(__expGenes);
        this.wildTypes = Utils.copy(__wildTypes);
        
    }
    
    public void setLatch(CountDownLatch latch) {
        this.latch = latch;
    }
    
    @Override
    public void run() {
        //synchronized(ComputeUnit.class) {
        InferBN.calScore(net, states, base_data, expGenes, wildTypes);
        //}
        
        latch.countDown();
    }
    
}
