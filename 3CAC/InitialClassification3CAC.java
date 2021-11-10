
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileWriter;
import java.util.ArrayList;
import java.util.Scanner;

public class InitialClassification3CAC {
	
	public static int binarySearch(ArrayList<Integer> list, int a){// return the index if exists, and return -(insertPos+1) if doesn't exist
		int l=0, r=list.size() - 1;
        while (l <= r) {
            int m = (l + r)/2;
            if (list.get(m) == a){
            	return m;
            }
            else if (list.get(m) < a){
            	l = m + 1;
            }
            else{
            	r = m - 1;
            }
        }
        return -1;
	}
	
	public static int binaryInsert(ArrayList<Integer> list, int a){
		int l=0, r=list.size() - 1;
        while (l <= r) {
            int m = (l + r)/2;
            if (list.get(m) == a){
            	return -1;
            }
            else if (list.get(m) < a){
            	l = m + 1;
            }
            else{
            	r = m - 1;
            }
        }
        return l;
	}
	
	public static int getContigID(String contig, boolean flyeAssembly){
		if(flyeAssembly){
			return Integer.parseInt(contig.substring(7));
		}
		else{
			String[] temp=contig.trim().split("_");
			return Integer.parseInt(temp[1]);
		}
	}
	
	public static void readPlasClass(String file_plasClass, ArrayList<Integer> contigs_plas, ArrayList<Double> score_plas, boolean flyeAssembly){
		try{
			Scanner in=new Scanner(new File(file_plasClass));
			String[] readLine; int contig_ID=-1, insert_pos=-1;
			
			while(in.hasNextLine()){
				readLine = in.nextLine().trim().split("[\\p{Space}]+");
				contig_ID = getContigID(readLine[0], flyeAssembly);
					
				insert_pos = binaryInsert(contigs_plas, contig_ID);
				contigs_plas.add(insert_pos, contig_ID);
				score_plas.add(insert_pos, Double.parseDouble(readLine[1]));
			}
			in.close();
			System.out.println(contigs_plas.size()+"  contigs are predicted by plasClass ! ");
			
		}catch(Exception e){
			
		}
	}
	
	public static void readDeepVirFinder(String file_deepVF, ArrayList<Integer> contigs_phage, ArrayList<Double> score_phage, ArrayList<Double> pvalue_phage, boolean flyeAssembly){
		try{
			Scanner in=new Scanner(new File(file_deepVF)); in.nextLine();
			String[] readLine; int contig_ID=-1, insert_pos=-1;
			
			while(in.hasNextLine()){
				readLine=in.nextLine().trim().split("[\\p{Space}]+");
				contig_ID = getContigID(readLine[0], flyeAssembly);
				
				insert_pos=binaryInsert(contigs_phage, contig_ID);
				contigs_phage.add(insert_pos, contig_ID);
				score_phage.add(insert_pos, Double.parseDouble(readLine[2]));
				pvalue_phage.add(insert_pos, Double.parseDouble(readLine[3]));
			}
			in.close();
			System.out.println(contigs_phage.size()+"  contigs are predicted by deepVirFinder ! ");
			
		}catch(Exception e){
			
		}
	}
	
	public static void readViralVerify(String file_outdir, String file_threeClass, ArrayList<Integer> contigs_plas, ArrayList<Double> score_plas, ArrayList<Integer> contigs_phage, ArrayList<Double> score_phage, ArrayList<Double> pvalue_phage, boolean flyeAssembly){
		try{
			BufferedWriter writer = new BufferedWriter(new FileWriter(file_outdir+"3CAC_initial_classification.fasta"));
			Scanner in=new Scanner(new File(file_threeClass)); 
			String[] readLine; int contig_ID=-1, insert_pos=-1; String contig_cate="";
			int num_phageC=0, num_plasC=0, num_uncertain=0;
			
			while(in.hasNextLine()) {
				readLine=in.nextLine().trim().split(",");
				contig_ID = getContigID(readLine[0], flyeAssembly);
				
				if(readLine[1].equalsIgnoreCase("Virus")) {//contigs classified as viruses by viralVerify
					insert_pos=binarySearch(contigs_phage, contig_ID);
					
					if(insert_pos!=-1){
						if(pvalue_phage.get(insert_pos)<=0.03){
							contig_cate="phage";
						}
						else if(score_phage.get(insert_pos)<=0.5){
							contig_cate="chromosome"; num_phageC++;
						}
						else{
							contig_cate="uncertain"; num_uncertain++;
						}
					}
					else{
						contig_cate="uncertain"; num_uncertain++;
					}
				}
				else if(readLine[1].equalsIgnoreCase("Plasmid")) {//contigs classified as plasmids by viralVerify
					insert_pos=binarySearch(contigs_plas, contig_ID);
					if(insert_pos!=-1){
						if(score_plas.get(insert_pos)>=0.7){
							contig_cate="plasmid";
						}
						else if(score_plas.get(insert_pos)<=0.3){
							contig_cate="chromosome"; num_plasC++;
						}
						else{
							contig_cate="uncertain"; num_uncertain++;
						}
					}
					else{
						System.out.println(readLine[0]+" classified as plasmid by viralVerify has no prediction of PlasClass !");
						contig_cate="uncertain"; num_uncertain++;
					}
				}
				else if(readLine[1].equalsIgnoreCase("Chromosome")) {//contigs classified as chromosomes by viralVerify
					contig_cate="chromosome";
				}
				else{
					contig_cate="uncertain";
				}
				writer.write(readLine[0]+","+contig_cate); writer.newLine();
			}
			in.close(); writer.close();
			
			System.out.println("Finish generating the initial classification of 3CAC !");
			System.out.println(num_phageC+" contigs classified as phages by viralVerify were reclassified as chromosome contigs by deepVirFinder. ");
			System.out.println(num_plasC+" contigs classified as plasmids by viralVerify were reclassified as chromosome contigs by plasClass. ");
			System.out.println(num_uncertain+" contigs classified as phages or plasmids by viralVerify were reclassified as uncertain. ");
			
		}catch(Exception e){
			
		}
	}
	
	public static void readPPRMeta(String file_outdir, String file_threeClass, ArrayList<Integer> contigs_plas, ArrayList<Double> score_plas, ArrayList<Integer> contigs_phage, ArrayList<Double> score_phage, ArrayList<Double> pvalue_phage, boolean flyeAssembly){
		try{
			BufferedWriter writer = new BufferedWriter(new FileWriter(file_outdir+"3CAC_initial_classification.fasta"));
			Scanner in=new Scanner(new File(file_threeClass)); in.nextLine();
			String[] readLine; int contig_ID=-1, insert_pos=-1; String contig_cate="";
			int num_phageC=0, num_plasC=0, num_uncertain=0;
			
			while(in.hasNextLine()) {
				readLine=in.nextLine().trim().split(",");
				contig_ID = getContigID(readLine[0], flyeAssembly);
				
				if(readLine[5].equalsIgnoreCase("phage") && Double.parseDouble(readLine[2])>=0.7) {//contigs classified as viruses by PPRMeta
					insert_pos=binarySearch(contigs_phage, contig_ID);
					
					if(insert_pos!=-1){
						if(pvalue_phage.get(insert_pos)<=0.03){
							contig_cate="phage";
						}
						else if(score_phage.get(insert_pos)<=0.5){
							contig_cate="chromosome"; num_phageC++;
						}
						else{
							contig_cate="uncertain"; num_uncertain++;
						}
					}
					else{
						contig_cate="uncertain"; num_uncertain++;
					}
				}
				else if(readLine[5].equalsIgnoreCase("plasmid") && Double.parseDouble(readLine[4])>=0.7) {//contigs classified as plasmids by PPRMeta
					insert_pos=binarySearch(contigs_plas, contig_ID);
					if(insert_pos!=-1){
						if(score_plas.get(insert_pos)>=0.7){
							contig_cate="plasmid";
						}
						else if(score_plas.get(insert_pos)<=0.3){
							contig_cate="chromosome"; num_plasC++;
						}
						else{
							contig_cate="uncertain"; num_uncertain++;
						}
					}
					else{
						System.out.println(readLine[0]+" classified as plasmid by PPRMeta has no prediction of PlasClass !");
						contig_cate="uncertain"; num_uncertain++;
					}
				}
				else if(readLine[5].equalsIgnoreCase("chromosome") && Double.parseDouble(readLine[3])>=0.7) {//contigs classified as chromosomes by PPRMeta
					contig_cate="chromosome";
				}
				else{
					contig_cate="uncertain";
				}
				writer.write(readLine[0]+","+contig_cate); writer.newLine();
			}
			in.close(); writer.close();
			
			System.out.println("Finish generating the initial classification of 3CAC !");
			System.out.println(num_phageC+" contigs classified as phages by PPRMeta were reclassified as chromosome contigs by deepVirFinder. ");
			System.out.println(num_plasC+" contigs classified as plasmids by PPRMeta were reclassified as chromosome contigs by plasClass. ");
			System.out.println(num_uncertain+" contigs classified as phages or plasmids by PPRMeta were reclassified as uncertain. ");
			
		}catch(Exception e){
			
		}
	}
	
	public static void main (String args[]){
		try{
			String assembler=""; //assembled by Flye or SPAdes
			String classifier=""; //use viralVerify or PPRMeta to generate the initial 3-way classification
			String file_threeClass=""; //three-class classification
			String file_plasClass=""; //plasClass classification
			String file_deepVF=""; //deepVirFinder classification
			String file_outdir=""; //output directory
			
			for(int i=0; i<10; i=i+2){
				if(args[i].equalsIgnoreCase("--assembler")){
					assembler=args[i+1];
					System.out.println("assembler: "+assembler);
				}
				else if(args[i].equalsIgnoreCase("--PPRMeta")){
					classifier="PPRMeta";
					file_threeClass=args[i+1];
					System.out.println("3-way classification file: "+file_threeClass);
				}
				else if(args[i].equalsIgnoreCase("--viralVerify")){
					classifier="viralVerify";
					file_threeClass=args[i+1];
					System.out.println("3-way classification file: "+file_threeClass);
				}
				else if(args[i].equalsIgnoreCase("--plasClass")){
					file_plasClass=args[i+1];
					System.out.println("PlasClass classification file: "+file_plasClass);
				}
				else if(args[i].equalsIgnoreCase("--deepVirFinder")){
					file_deepVF=args[i+1];
					System.out.println("deepVirFiner classification file: "+file_deepVF);
				}
				else if(args[i].equalsIgnoreCase("--output")){
					file_outdir=args[i+1];
					if(file_outdir.charAt(file_outdir.length()-1) != '/'){
						file_outdir=file_outdir+"/";
					}
					System.out.println("output directory: "+file_outdir);
				}
			}
			System.out.println();
			
			boolean flyeAssembly=false;  boolean viralVerify=false;
			if(assembler.equalsIgnoreCase("flye") || assembler.equalsIgnoreCase("metaflye")){
				flyeAssembly=true;
			}
			if(classifier.equalsIgnoreCase("viralVerify")){
				viralVerify=true;
			}
			
			//read plasClass classification
			ArrayList<Integer> contigs_plas=new ArrayList<Integer>();
			ArrayList<Double> score_plas=new ArrayList<Double>();
			readPlasClass(file_plasClass, contigs_plas, score_plas, flyeAssembly);
			
			
			//read deepVirFinder classification
			ArrayList<Integer> contigs_phage=new ArrayList<Integer>();
			ArrayList<Double> score_phage=new ArrayList<Double>();
			ArrayList<Double> pvalue_phage=new ArrayList<Double>();
			readDeepVirFinder(file_deepVF, contigs_phage, score_phage, pvalue_phage, flyeAssembly);
			
			//read three-way classification
			if(viralVerify){
				readViralVerify(file_outdir, file_threeClass, contigs_plas, score_plas, contigs_phage, score_phage, pvalue_phage, flyeAssembly);
			}
			else{
				readPPRMeta(file_outdir, file_threeClass, contigs_plas, score_plas, contigs_phage, score_phage, pvalue_phage, flyeAssembly);
			}
			System.out.println("The initial classification of 3CAC can be found in: "+file_outdir+"3CAC_initial_classification.fasta");
			
		}catch(Exception e){
			
		}
	}

}

