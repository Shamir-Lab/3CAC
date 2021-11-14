
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileWriter;
import java.util.ArrayList;
import java.util.Scanner;

public class Classify3CAC {
	
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
	
	public static int getContigID(String contig, boolean flyeAssembly){
		if(flyeAssembly){
			if(contig.substring(0, 6).equalsIgnoreCase("contig")){
				return Integer.parseInt(contig.substring(7));
			}
			else{
				return Integer.parseInt(contig.substring(9));
			}
		}
		else{
			String[] temp=contig.trim().split("_");
			return Integer.parseInt(temp[1]);
		}
	}
	
	public static int getEdgeID(String edge, boolean flyeAssembly){
		if(flyeAssembly){
			return Integer.parseInt(edge.substring(5));
		}
		else{
			return Integer.parseInt(edge);
		}
	}
	
	public static int flyeContigsToEdgesInGraph(String file_assembly_info, ArrayList<Integer> contigs, ArrayList<ArrayList<Integer>> contigEdges) {
		try{
			//read edges of each contig
			int max_edge_ID=0;
			Scanner in=new Scanner(new File(file_assembly_info)); in.nextLine();  
			String[] readLine, edgeInPath;
			int contig_ID=-1, insert_pos=-1, edge_ID=-1;
			while(in.hasNextLine()) {
				readLine=in.nextLine().trim().split("[\\p{Space}]+");
				edgeInPath=readLine[7].trim().split(",");
				ArrayList<Integer> edgesInCurContig=new ArrayList<Integer>();
				for(int i=0;i<edgeInPath.length;i++) {
					if(!edgeInPath[i].equalsIgnoreCase("*") && !edgeInPath[i].equalsIgnoreCase("??")) {
						edge_ID = Math.abs(Integer.parseInt(edgeInPath[i]));
						edgesInCurContig.add(edge_ID);
						if(max_edge_ID < edge_ID){
							max_edge_ID = edge_ID;
						}
					}
				} 
				if(readLine[0].substring(0, 6).equalsIgnoreCase("contig")){
					contig_ID=Integer.parseInt(readLine[0].substring(7));
				}
				else{
					contig_ID=Integer.parseInt(readLine[0].substring(9));
				}
				insert_pos=binaryInsert(contigs, contig_ID);
				
				contigs.add(insert_pos, contig_ID);
				contigEdges.add(insert_pos, edgesInCurContig);
			}
			in.close();
			System.out.println(contigs.size()+"  contigs contained in "+file_assembly_info+" !");
			return max_edge_ID;
			
		}catch(Exception e){
			
		}
		return -1;
	}	
	
	public static int spadesContigsToEdgesInGraph(String file_assembly_info, ArrayList<Integer> contigs, ArrayList<ArrayList<Integer>> contigEdges) {
		try{
			//read edges of each contig
			int max_edge_ID=0;
			Scanner in=new Scanner(new File(file_assembly_info));   
			ArrayList<Integer> edgesInCurContig=new ArrayList<Integer>();
			int contig_ID=-1, edge_ID=-1, insert_pos=-1, length=-1; boolean readNode=false;
			String readLine; String[] nodeLabel, edgeInPath;
			
			while(in.hasNextLine()) {
				readLine=in.nextLine();
				if(readLine.length()>4 && readLine.substring(0, 4).equalsIgnoreCase("NODE")){//read contig ID
					length=readLine.length();
					if(readLine.charAt(length-1)=='\''){
						readNode=false;
					}
					else{
						readNode=true;
						if(edgesInCurContig.size()>0){
							ArrayList<Integer> temp=new ArrayList<Integer>();
							for(int i=0;i<edgesInCurContig.size();i++){
								temp.add(edgesInCurContig.get(i));
							}
							contigEdges.add(insert_pos, temp);
							edgesInCurContig.clear();
						}
						
						nodeLabel=readLine.trim().split("_");
						contig_ID=Integer.parseInt(nodeLabel[1]);
						insert_pos=binaryInsert(contigs, contig_ID);
						contigs.add(insert_pos, contig_ID);
						
					}
				}
				else{
					if(readNode){//read the path if this Node was add to contigs
						length=readLine.length();
						if(readLine.charAt(length-1)==';'){
							readLine=readLine.substring(0, length-1);
						}
						edgeInPath=readLine.trim().split(",");
						for(int i=0;i<edgeInPath.length;i++){
							length=edgeInPath[i].length();
							edge_ID=Integer.parseInt(edgeInPath[i].substring(0, length-1));
							edgesInCurContig.add(edge_ID);
							
							if(max_edge_ID < edge_ID){
								max_edge_ID = edge_ID;
							}
						}
					}
				}
			}
			in.close();
			
			if(edgesInCurContig.size()>0){
				contigEdges.add(insert_pos, edgesInCurContig);
			}
			System.out.println(contigs.size()+"  contigs contained in "+file_assembly_info+" !");
			return max_edge_ID;
			
		}catch(Exception e){
			
		}
		return -1;
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
	
	
	public static void initialClassification(String file_outdir, boolean viralVerify, boolean flyeAssembly, String file_threeClass, ArrayList<Integer> contigs_plas, ArrayList<Double> score_plas, ArrayList<Integer> contigs_phage, ArrayList<Double> score_phage, ArrayList<Double> pvalue_phage, ArrayList<String> initial_class){
		try{
			BufferedWriter writer = new BufferedWriter(new FileWriter(file_outdir+"3CAC_initial_classification.fasta"));
			Scanner in=new Scanner(new File(file_threeClass)); 
			if(!viralVerify){//ignore the header line for PPR-Meta
				in.nextLine();
			}
			
			String[] readLine; int contig_ID=-1, contig_cate=-1, insert_pos=-1; 
			int num_phageC=0, num_plasC=0, num_uncertain=0;
			
			while(in.hasNextLine()) {
				readLine=in.nextLine().trim().split(",");
				contig_ID = getContigID(readLine[0], flyeAssembly);
				contig_cate=getCate(readLine, viralVerify);
				
				if(contig_cate == 1) {//contigs classified as viruses by viralVerify (or PPRMeta)
					insert_pos=binarySearch(contigs_phage, contig_ID);
					if(insert_pos >= 0){
						if(pvalue_phage.get(insert_pos)<=0.03){
							contig_cate=1;
						}
						else if(score_phage.get(insert_pos)<=0.5){
							contig_cate=3; num_phageC++;
						}
						else{
							contig_cate=4; num_uncertain++;
						}
					}
					else{
						contig_cate=4; num_uncertain++;
					}
				}
				else if(contig_cate == 2) {//contigs classified as plasmids by viralVerify (or PPRMeta)
					insert_pos=binarySearch(contigs_plas, contig_ID);
					if(insert_pos >= 0){
						if(score_plas.get(insert_pos)>=0.7){
							contig_cate=2;
						}
						else if(score_plas.get(insert_pos)<=0.3){
							contig_cate=3; num_plasC++;
						}
						else{
							contig_cate=4; num_uncertain++;
						}
					}
					else{
						System.out.println(readLine[0]+" classified as plasmid by viralVerify has no prediction of PlasClass !");
						contig_cate=4; num_uncertain++;
					}
				}
				
				initial_class.add(readLine[0]+","+contig_cate);
				writer.write(readLine[0]+","+contig_cate); writer.newLine();
			}
			in.close(); writer.close();
			contigs_plas.clear(); score_plas.clear(); contigs_phage.clear(); score_phage.clear(); pvalue_phage.clear();
			System.out.println("Finish generating the initial classification of 3CAC !");
			System.out.println(num_phageC+" contigs classified as phages by viralVerify were reclassified as chromosome contigs by deepVirFinder. ");
			System.out.println(num_plasC+" contigs classified as plasmids by viralVerify were reclassified as chromosome contigs by plasClass. ");
			System.out.println(num_uncertain+" contigs classified as phages or plasmids by viralVerify were reclassified as uncertain. ");
			
		}catch(Exception e){
			
		}
	}
		
	public static void contigCategoryToEdges(ArrayList<String> initial_class, ArrayList<Integer> contigs, ArrayList<ArrayList<Integer>> contigEdges, boolean flyeAssembly, int[] edgesCate){		
		
		String[] readLine; int contig_ID=-1, contig_cate=-1, contig_pos=-1, edge_ID=-1, edge_insert_pos=-1;
		for(int i=0; i<initial_class.size(); i++) {
			readLine = initial_class.get(i).trim().split(",");
			contig_cate = Integer.parseInt(readLine[1]);
			
			//assign category of contigs to edges
			if(contig_cate !=4){
				contig_ID = getContigID(readLine[0], flyeAssembly);
				contig_pos = binarySearch(contigs, contig_ID);
				for(int j=0; j<contigEdges.get(contig_pos).size();j++) {
					edge_ID=contigEdges.get(contig_pos).get(j)-1;
					
					if(edgesCate[edge_ID] == 4){
						edgesCate[edge_ID]=contig_cate;
					}
					else if(edgesCate[edge_ID] != contig_cate){
						edgesCate[edge_ID]=-1; //edge with conflict categories as it was contained in multiple contigs with different categories
					}
				}
			}
		}
		
		for(int i=0; i<edgesCate.length; i++){
			if(edgesCate[i]==-1){ //set conflict edges to be uncertain
				edgesCate[i]=4; 
			}
		}
		System.out.println("Edges contained in contigs are assigned with contig categories ! ");
		
	}
	
	public static void readAssemblyGraph(String file_assembly_graph, boolean flyeAssembly, int[] edgesLength, ArrayList<ArrayList<Integer>> edgesAdj) {
		try{
			System.out.println("Reading assembly graph "+file_assembly_graph+" !");
			
			Scanner in=new Scanner(new File(file_assembly_graph));   
			if(flyeAssembly){
				in.nextLine();
			}
			
			String[] readLine; boolean end=false; 
			int edge_ID=-1, insert_pos=-1, preEdge_ID=-1, sufEdge_ID=-1;
			while(in.hasNextLine() && !end) {
				readLine=in.nextLine().trim().split("[\\p{Space}]+");
				if(readLine[0].equalsIgnoreCase("S")) {
					edge_ID=getEdgeID(readLine[1], flyeAssembly)-1;
					edgesLength[edge_ID] = readLine[2].length();
				}
				else if(readLine[0].equalsIgnoreCase("L")) {
					if(!readLine[1].equalsIgnoreCase(readLine[3])){
						preEdge_ID = getEdgeID(readLine[1], flyeAssembly)-1;
						sufEdge_ID = getEdgeID(readLine[3], flyeAssembly)-1;
						
						insert_pos=binarySearch(edgesAdj.get(preEdge_ID), sufEdge_ID);
						if(insert_pos < 0){
							insert_pos = Math.abs(insert_pos)-1;
							edgesAdj.get(preEdge_ID).add(insert_pos, sufEdge_ID);
							
							insert_pos=binarySearch(edgesAdj.get(sufEdge_ID), preEdge_ID);
							insert_pos = Math.abs(insert_pos)-1;
							edgesAdj.get(sufEdge_ID).add(insert_pos, preEdge_ID);
						}
					}
				}
				else {
					end=true;
				}
			}
			in.close();
			System.out.println("End Reading assembly graph !");
			
		}catch(Exception e){
			
		}
	}

	public static void correction(ArrayList<ArrayList<Integer>> edgesAdj, int[] edgesCate){
		
		int curEdgeCate=-1, curEdgeAdj=-1, adjCate=-1, correctEdges=0; 
		int num_support=0, num_same=0, loop_num=0;  
		boolean multiAdjCates=false, correct=true; 
		
		while(correct){
			correct=false; loop_num++;
			for(int i=0;i<edgesAdj.size();i++){
				curEdgeCate=edgesCate[i];
				if(curEdgeCate!=4 && edgesAdj.get(i).size() >= 2){// get a classified edge with at least two neighbors
					
					num_same=0; multiAdjCates=false; 
					for(int j=0; j<edgesAdj.get(i).size() && !multiAdjCates; j++){//check adjacent edges for current edge
						curEdgeAdj=edgesAdj.get(i).get(j);
						adjCate=edgesCate[curEdgeAdj];
						
						if (adjCate == curEdgeCate){
							num_same++;
						}
						else if(adjCate != 4){//find a classified neighboring edge with different category
							num_support=1; 
							j++;
							
							while(!multiAdjCates && j<edgesAdj.get(i).size()){
								curEdgeAdj=edgesAdj.get(i).get(j);
								if(edgesCate[curEdgeAdj] != 4){
									if(adjCate != edgesCate[curEdgeAdj]){
										multiAdjCates=true;
									}
									else{
										num_support++;
									}
								}
								j++;
							}
							
							//correct current edge category to be same as its classified neighbors
							if(!multiAdjCates && num_same==0 && num_support>=2){
								edgesCate[i] = adjCate;
								correct=true;
								correctEdges++;
							}
						}
					}
					
				}
			}
		}
		System.out.println("The correction step was repeated for "+loop_num+" times, and "+correctEdges+"  classified edges were corrected !");
//		System.out.println("finish correction steps !  ");
	}
	
    public static void propagation(ArrayList<ArrayList<Integer>> edgesAdj, int[] edgesCate){
		
		int curEdgeCate=-1, curEdgeAdj=-1, adjCate=-1, propagateEdges=0; 
		int num_support=0, num_same=0, loop_num=0;  
		boolean multiAdjCates=false, propagate=true; 
		
		while(propagate){
			propagate=false; loop_num++;
			for(int i=0;i<edgesAdj.size();i++){
				if(edgesCate[i]==4 && edgesAdj.get(i).size() >= 1){ //get an uncertain edge with at least one neighbor
					
					multiAdjCates=false; adjCate=-1; 
					for(int j=0; j<edgesAdj.get(i).size() && !multiAdjCates; j++){// check neighbors for current edge
						curEdgeAdj=edgesAdj.get(i).get(j);
						if(edgesCate[curEdgeAdj] != 4){//find a classified neighbors
							if(adjCate==-1){
								adjCate=edgesCate[curEdgeAdj]; 
							}
							else{
								if(adjCate != edgesCate[curEdgeAdj]){
									multiAdjCates=true;
								}
							}
						}												
					}
					
					if(!multiAdjCates && adjCate!=-1){ // propagate category of classified neighbors to the uncertain edge
						edgesCate[i] = adjCate;
						propagate=true;
						propagateEdges++;
					}
				}
			}
		}
		
		System.out.println("The propagation step was repeated for "+loop_num+" times, and "+propagateEdges+"  uncertain edges were classified !");
//		System.out.println("finish propagation steps !  ");
	}
    
    public static void finalContigCategoryFromEdges(String file_outdir, ArrayList<Integer> contigs, ArrayList<ArrayList<Integer>> contigEdges, int[] edgesLength, int[] edgesCate, boolean flyeAssembly){
    	try{
    		BufferedWriter writer = new BufferedWriter(new FileWriter(file_outdir+"3CAC_classification.fasta"));
    		int edge_ID=-1, edge_cate=-1, final_cate=-1; double sumLength=0, maxLength=0;
    		int[] cate_length=new int[4];
            String prefix="NODE_"; String finalContigCate="";
    		if(flyeAssembly){
    			prefix="contig_";
    		}
    		
    		for(int i=0; i<contigs.size(); i++){
    			for(int k=0;k<cate_length.length;k++){
    				cate_length[k]=0;
    			}
    			
    			maxLength=0; sumLength=0;
    			for(int j=0; j<contigEdges.get(i).size(); j++){
    				edge_ID=contigEdges.get(i).get(j)-1;
    				edge_cate=edgesCate[edge_ID]-1;
    				cate_length[edge_cate]=cate_length[edge_cate]+edgesLength[edge_ID];
    				sumLength=sumLength+edgesLength[edge_ID];
    				
    				if(maxLength<cate_length[edge_cate]){
    					maxLength=cate_length[edge_cate]; final_cate=edge_cate+1;
    				}
    			}
    			
    			//for contigs contain edges from different categories.
    			//if the categories with longest length doesn't account for 80% of the contig length, this contig will be classified as uncertain
    			if(maxLength/sumLength<0.8){ 
    				final_cate=4;
    			}
    			
    			if(final_cate==1){
    				finalContigCate="phage";
    			}
    			else if(final_cate==2){
    				finalContigCate="plasmid";
    			}
    			else if(final_cate==3){
    				finalContigCate="chromosome";
    			}
    			else{
    				finalContigCate="uncertain";
    			}
    			writer.write(prefix+contigs.get(i)+","+finalContigCate); writer.newLine();
    		}
    		writer.close();
    		System.out.println("The classification result of 3CAC can be found in: "+file_outdir+"3CAC_classification.fasta");
    		
    	}catch(Exception e){
			
		}
    }
		
	public static void main (String args[]){
		try{
			String assembler=""; // assembled by Flye or SPAdes
			String file_outdir=""; // output directory
			String file_path=""; // file containing path information for each contig
			String file_graph=""; // assembly graph file
			String file_threeClass=""; //three-class classification (viralVerify or PPRMeta)
			String file_plasClass=""; //plasClass classification
			String file_deepVF=""; //deepVirFinder classification
			boolean viralVerify=false;
			
			for(int i=0; i<14; i=i+2){
				if(args[i].equalsIgnoreCase("--assembler")){
					assembler=args[i+1];
				}
				else if(args[i].equalsIgnoreCase("--viralVerify")){
					viralVerify=true;
					file_threeClass=args[i+1];
				}
				else if(args[i].equalsIgnoreCase("--PPRMeta")){
					viralVerify=false;
					file_threeClass=args[i+1];
				}
				else if(args[i].equalsIgnoreCase("--plasClass")){
					file_plasClass=args[i+1];
				}
				else if(args[i].equalsIgnoreCase("--deepVirFinder")){
					file_deepVF=args[i+1];
				}
				else if(args[i].equalsIgnoreCase("--graph")){
					file_graph=args[i+1];
				}
				else if(args[i].equalsIgnoreCase("--path")){
					file_path=args[i+1];
				}
				else if(args[i].equalsIgnoreCase("--output")){
					file_outdir=args[i+1];
					if(file_outdir.charAt(file_outdir.length()-1) != '/'){
						file_outdir=file_outdir+"/";
					}
				}
			}

			System.out.println();
			System.out.println("Assembler: "+assembler);
			if(viralVerify){
				System.out.println("viralVerify classification: "+file_threeClass);
			}
			else{
				System.out.println("PPR-Meta classification: "+file_threeClass);
			}
			System.out.println("PlasClass classification: "+file_plasClass);
			System.out.println("DeepVirFinder classification: "+file_deepVF);
			System.out.println("assembly graph file: "+file_graph);
			System.out.println("contig path file: "+file_path);
			System.out.println("output directory: "+file_outdir);
			System.out.println();
			
			boolean flyeAssembly=false; 
			if(assembler.equalsIgnoreCase("flye") || assembler.equalsIgnoreCase("metaflye")){
				flyeAssembly=true;
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
			
			
			//generate initial classification of 3CAC
			ArrayList<String> initialClass=new ArrayList<String>();
			initialClassification(file_outdir, viralVerify, flyeAssembly, file_threeClass, contigs_plas, score_plas, contigs_phage, score_phage, pvalue_phage, initialClass);
			
			
			//read edges in contigs
			int max_edge_ID=0;
			ArrayList<Integer> contigs=new ArrayList<Integer>();
			ArrayList<ArrayList<Integer>> contigEdges=new ArrayList<ArrayList<Integer>>();
			if(flyeAssembly){
				max_edge_ID = flyeContigsToEdgesInGraph(file_path, contigs, contigEdges);
			}
			else{
				max_edge_ID = spadesContigsToEdgesInGraph(file_path, contigs, contigEdges);
			}
			
			
			//read assembly graph to get adjacency and edge length
			int[] edgesLength=new int[max_edge_ID];
			int[] edgesCate=new int[max_edge_ID];
			ArrayList<ArrayList<Integer>> edgesAdj=new ArrayList<ArrayList<Integer>>();
			for(int i=0; i<edgesLength.length; i++){
				edgesLength[i]=-1; edgesCate[i]=4;
				ArrayList<Integer> tempAdj=new ArrayList<Integer>();
				edgesAdj.add(tempAdj);
			}
			readAssemblyGraph(file_graph, flyeAssembly, edgesLength, edgesAdj);
			
			
			//assign contig categories to edges 
			contigCategoryToEdges(initialClass, contigs, contigEdges, flyeAssembly, edgesCate);
			
			
			//correction and propagation step
			correction(edgesAdj, edgesCate);
			propagation(edgesAdj, edgesCate);
			for(int i=0; i<edgesAdj.size(); i++){
				edgesAdj.get(i).clear();
			}
			edgesAdj.clear();
			
			
			//final contig category decided by edge category
			finalContigCategoryFromEdges(file_outdir, contigs, contigEdges, edgesLength, edgesCate, flyeAssembly);
			
			System.out.println("Finish 3CAC algorithm !");
    		
		}catch(Exception e){
			
		}
	}

}

