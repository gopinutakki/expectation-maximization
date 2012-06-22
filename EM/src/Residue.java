/**
 * Class: Residue
 * 
 * Definition for holding the residue counts, frequencies.
 * @author gopi
 *
 */
public class Residue {

	double dna[][] = null;
	double aa[][] = null;

	public Residue(int motifLen) {
		dna = new double[4][motifLen + 1];
		aa = new double[20][motifLen + 1];		
		resetDNA(motifLen);
	}
	
	public void resetDNA(int motifLen){
		for(int i = 0; i < 4; i++){
			for(int j = 0; j < motifLen+1; j++)
				dna[i][j] = 0;
		}
	}
}
