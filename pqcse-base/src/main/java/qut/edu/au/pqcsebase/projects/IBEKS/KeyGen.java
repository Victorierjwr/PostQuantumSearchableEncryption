package qut.edu.au.pqcsebase.projects.IBEKS;

import qut.edu.au.pqcsebase.projects.IBEKS.vo.KeyPair;
import qut.edu.au.pqcsebase.projects.IBEKS.vo.MSK;
import qut.edu.au.pqcsebase.tools.BigIntMatrix;


import static qut.edu.au.pqcsebase.algorithms.BasisDel.basisDel;
import static qut.edu.au.pqcsebase.projects.IBEKS.Params.n;
import static qut.edu.au.pqcsebase.projects.IBEKS.Params.sigma;
import static qut.edu.au.pqcsebase.tools.BigIntMatrix.inverseQ;

public class KeyGen {

    public static KeyPair keyGen(Setup setup, MSK msk, String identity) {

        BigIntMatrix A = setup.pp.A;
        BigIntMatrix TA = msk.msk;

        BigIntMatrix R_id = setup.hash1(identity);
        BigIntMatrix A_id = A.multiplyQ(inverseQ(R_id));

        BigIntMatrix T_id = basisDel(A,R_id,TA,sigma,n);

        return new KeyPair(A_id,T_id);
    }
}
