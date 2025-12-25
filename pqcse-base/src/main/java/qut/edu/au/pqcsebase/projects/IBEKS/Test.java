package qut.edu.au.pqcsebase.projects.IBEKS;

import qut.edu.au.pqcsebase.projects.IBEKS.vo.CT;
import qut.edu.au.pqcsebase.projects.IBEKS.vo.TW;

import java.math.BigInteger;

import static qut.edu.au.pqcsebase.projects.IBEKS.Params.q;
import static qut.edu.au.pqcsebase.tools.VectorUtils.innerProduct;

public class Test {

    public static Boolean test(CT ct, TW tw) {
        BigInteger e0cw = innerProduct(tw.getE0(), ct.getC_w());
        BigInteger e1cid_e0cw = BigInteger.ZERO;
        int len = tw.getE1().length;
        for (int i = 0; i < tw.getE1().length; i++) {
            if(i < len - 1){
                e1cid_e0cw = e1cid_e0cw.add(tw.getE1()[i].multiply(ct.getC_id()[i]));
            }else{
                e1cid_e0cw = e1cid_e0cw.add(tw.getE1()[i].multiply(e0cw));
            }
        }
        System.out.println("e1cid_e0cw: " + e1cid_e0cw);
        System.out.println("e1cid_e0cw mod q : " + e1cid_e0cw.mod(BigInteger.valueOf(q)));
        BigInteger u = (ct.getC_u().subtract(e1cid_e0cw)).mod(BigInteger.valueOf(q));
        System.out.println("u: " + u);
        return u.intValue() <= Math.ceil((double) q / 4);
    }
}
