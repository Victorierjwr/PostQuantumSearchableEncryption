package qut.edu.au.pqcsebaseopt.projects.IBEKS;



import qut.edu.au.pqcsebaseopt.projects.IBEKS.vo.MSK;
import qut.edu.au.pqcsebaseopt.projects.IBEKS.vo.PP;
import qut.edu.au.pqcsebaseopt.tools.BigIntMatrix;
import qut.edu.au.pqcsebaseopt.tools.LatticeHashUtils;

import java.math.BigInteger;
import java.util.List;

import static qut.edu.au.pqcsebaseopt.algorithms.TrapGen.trapGen;
import static qut.edu.au.pqcsebaseopt.projects.IBEKS.Params.*;
import static qut.edu.au.pqcsebaseopt.tools.VectorUtils.innerProduct;
import static qut.edu.au.pqcsebaseopt.tools.VectorUtils.randomVector;


public class Setup {

    public PP pp;
    public MSK msk;
    
    public long m;


    private LatticeHashUtils hashUtils;

    public Setup(){
        BigIntMatrix[] matrices = trapGen(n,q,flag);
        this.pp = new PP(); this.msk = new MSK(matrices[1]);
        this.pp.A = matrices[0]; m = this.pp.A.getCols();
        this.pp.B = BigIntMatrix.random((int)n,(int)m, BigInteger.valueOf(q));
        this.pp.u = randomVector((int)(n), BigInteger.valueOf(q));
        List<BigInteger[]> innerProd = randInnerProduct();
        this.pp.b1 = innerProd.get(0);
        this.pp.y1 = innerProd.get(1);
        this.hashUtils = new LatticeHashUtils(BigInteger.valueOf(q),(int)m);
    }

    public BigIntMatrix hash1(String str){
        return this.hashUtils.StringHashToMatrixModQ(str);
    }

    public BigInteger hash2(String str){
        return this.hashUtils.StringHashToValueModQ(str);
    }

    private List<BigInteger[]> randInnerProduct(){
        BigInteger[] b1 = randomVector((int)(n_N), BigInteger.valueOf(q));
        BigInteger[] y1 = randomVector((int)(n_N), BigInteger.valueOf(q));
        if(innerProduct(b1, y1).mod(BigInteger.valueOf(q)).equals(BigInteger.ZERO)){
            return List.of(b1, y1);
        }else{
            return randInnerProduct();
        }
    }
}
