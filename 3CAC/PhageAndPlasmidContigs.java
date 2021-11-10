
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileWriter;
import java.util.ArrayList;
import java.util.Scanner;

public class PhageAndPlasmidContigs {
	
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
        return -(l+1);
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
	
	public static int getContigID(String header){	
		
		String[] index;
		if(header.substring(0, 4).equalsIgnoreCase("NODE")){// contigs assembled from SPAdes
			index=header.trim().split("_");
			return Integer.parseInt(index[1]);
		}
		else{ //contigs assembled from Flye
			return Integer.parseInt(header.substring(7));
		}
	}
	
	public static void readContigs(String file_contig, ArrayList<Integer> contigs_ID, ArrayList<ArrayList<String>> contigs_Seq){
		try{
			
			Scanner in=new Scanner(new File(file_contig)); 
			ArrayList<String> seqs=new ArrayList<String>();
			String readLine="";  int curContig_ID=-1, insert_pos=-1;
			 
			while(in.hasNextLine()) {
				readLine=in.nextLine();
				if(readLine.length()>0){
					if(readLine.charAt(0)=='>'){
						if(seqs.size()>0){
							ArrayList<String> curContig_Seq=new ArrayList<String>();
							for(int i=0; i<seqs.size(); i++){
								curContig_Seq.add(seqs.get(i));
							}
							contigs_Seq.add(insert_pos, curContig_Seq);
						}
						seqs.clear();
						
						curContig_ID=getContigID(readLine.substring(1));
						insert_pos=binaryInsert(contigs_ID, curContig_ID);
						
						contigs_ID.add(insert_pos, curContig_ID);
						seqs.add(readLine);
					}
					else{
						seqs.add(readLine);
					}
				}
			}
			in.close(); 
			if(seqs.size()>0){
				contigs_Seq.add(insert_pos, seqs);
			}
			System.out.println(contigs_ID.size()+" contigs are contained in file: "+file_contig);
			
			
		}catch(Exception e){
			
		}
	}
	
	public static int getCate(String[] readLine, boolean viralVerify){
		
		double PPRMeta_score=0.7;
		if(viralVerify){
			if(readLine[1].equalsIgnoreCase("Virus")){
				return 1;
			}
			else if(readLine[1].equalsIgnoreCase("Plasmid")){
				return 2;
			}
			else if(readLine[1].equalsIgnoreCase("Chromosome")){
				return 3;
			}
			else{
				return 4;
			}
		}
		else{
			if(readLine[5].equalsIgnoreCase("phage") && Double.parseDouble(readLine[2])>=PPRMeta_score){
				return 1;
			}
			else if(readLine[5].equalsIgnoreCase("plasmid") && Double.parseDouble(readLine[4])>=PPRMeta_score){
				return 2;
			}
			else if(readLine[5].equalsIgnoreCase("chromosome") && Double.parseDouble(readLine[3])>=PPRMeta_score){
				return 3;
			}
			else{
				return 4;
			}
		}
	}
	
    
	public static void writePhagesAndPlasmids(boolean viralVerify, String file_classification, ArrayList<Integer> contigs_ID, ArrayList<ArrayList<String>> contigs_Seq){
		try{
			BufferedWriter writer_phage = new BufferedWriter(new FileWriter("phageContigs.fasta"));
			BufferedWriter writer_plasmid = new BufferedWriter(new FileWriter("plasmidContigs.fasta"));
			
			Scanner in=new Scanner(new File(file_classification));
			if(!viralVerify){//ignore the header line for 
				in.nextLine();
			}
			
			String[] readLine; int contig_cate=-1, contig_ID=-1, insert_pos=-1; 
			int num_phage=0, num_plasmid=0;
			while(in.hasNextLine()){
				readLine=in.nextLine().trim().split(",");
				contig_cate=getCate(readLine, viralVerify);
				
				if(contig_cate==1){
					contig_ID=getContigID(readLine[0]);
					insert_pos=binarySearch(contigs_ID, contig_ID);
					
					for(int i=0; i<contigs_Seq.get(insert_pos).size(); i++){
						writer_phage.write(contigs_Seq.get(insert_pos).get(i)); writer_phage.newLine();
					}
					num_phage++;
				}
				else if(contig_cate==2){
					contig_ID=getContigID(readLine[0]);
					insert_pos=binarySearch(contigs_ID, contig_ID);
					
					for(int i=0; i<contigs_Seq.get(insert_pos).size(); i++){
						writer_plasmid.write(contigs_Seq.get(insert_pos).get(i)); writer_plasmid.newLine();
					}
					num_plasmid++;
				}
			}
			in.close(); writer_phage.close(); writer_plasmid.close();
			System.out.println(num_phage+" phage contigs and "+num_plasmid+" plasmid contigs were predicted in "+file_classification);
			
		}catch(Exception e){
			
		}
		
	}
	
	
	public static void main (String args[]){
		try{
			
			String file_contig=""; // file of contigs to be classified
			String file_classification=""; //classification result of PPR-Meta or viralVerify
			boolean viralVerify=false;
			
			for(int i=0; i<4; i=i+2){
				if(args[i].equalsIgnoreCase("--contig")){
					file_contig=args[i+1];
					System.out.println("contigs to be classified: "+file_contig);
				}
				else if(args[i].equalsIgnoreCase("--PPRMeta")){
					viralVerify=false;
					file_classification=args[i+1];
					System.out.println("PPRMeta classification file: "+file_classification);
				}
				else if(args[i].equalsIgnoreCase("--viralVerify")){
					viralVerify=true;
					file_classification=args[i+1];
					System.out.println("viralVerify classification file: "+file_classification);
				}
			}
			System.out.println();
			
			
			//read contig file
			ArrayList<Integer> contigs_ID=new ArrayList<Integer>();
			ArrayList<ArrayList<String>> contigs_Seq=new ArrayList<ArrayList<String>>();
			readContigs(file_contig, contigs_ID, contigs_Seq);
			
			
			//write phage contigs and plasmid contigs
			writePhagesAndPlasmids(viralVerify, file_classification, contigs_ID, contigs_Seq);
			
			
			
		}catch(Exception e){
			
		}
	}

}

