package bni.comp;
import java.util.ArrayList;
import java.util.concurrent.atomic.AtomicInteger;
import mod.Mutation;
import mod.jmut.core.comp.Node;

public class EdgeAddition2 implements Mutation {
    int pos;
    int type;
    
    public EdgeAddition2(int pos, int type) {
        this.pos = pos;
        this.type = type;
    }
    
    public void nodeMutation(ArrayList<Integer> input, ArrayList<Integer> I, 
                             ArrayList<Integer> O, AtomicInteger O_default) {
    }

    public void edgeMutation(int iSrc, ArrayList<Integer> tar_input, ArrayList<Integer> tar_I,
                             ArrayList<Integer> tar_O, AtomicInteger tar_O_default) {
        int Ik, Ok;
        int k = tar_input.size() + 1;
        double temp = Math.exp(0 - Node.THETA * (Math.pow(2, -k)));
        double probO = temp / (1 + temp);
        
//        Ik = (Math.random() < 0.5) ? 1 : 0;        
        Ok = (Math.random() < probO) ? 1 : 0;
        Ik = Ok;
        if(this.type == -1) Ik = 1 - Ok;
        
        //Insert the new edge into the first position in the input nodes
        tar_input.add(pos, iSrc);
        tar_I.add(pos, Ik);
        tar_O.add(pos, Ok);
        
        if(k == 1) {
            //Initially, the target node has no incoming links
            tar_O_default.set(1 - Ok);
        }
    }    
}