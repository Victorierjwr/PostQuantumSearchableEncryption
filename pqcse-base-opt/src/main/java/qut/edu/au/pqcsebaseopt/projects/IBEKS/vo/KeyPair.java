package qut.edu.au.pqcsebaseopt.projects.IBEKS.vo;

import lombok.*;
import qut.edu.au.pqcsebaseopt.tools.BigIntMatrix;

@Data
@Setter @Getter
@AllArgsConstructor@NoArgsConstructor
@ToString
public class KeyPair {

    public BigIntMatrix pk_id;
    public BigIntMatrix sk_id;

}
