#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <math.h>
#include <pthread.h>
#include <windows.h>
#include <assert.h>

int getNumberOfCoresWindows() {
    SYSTEM_INFO sysinfo;
    GetSystemInfo(&sysinfo);
    return sysinfo.dwNumberOfProcessors;
}

/* Use this if running on linux. Remember to comment out the other core detecting method.
int getNumberOfCores(){
	return sysconf(_SC_NPROCESSORS_ONLN);
}
*/

//------------------------------------------------------------------------------------------------------------------------------------------
struct kmer_s{
	char kmerChars[100]; 
	
};

struct node_s {
	struct kmer_s data;
	struct node_s *next;
};

struct hashtable_s {
	int size;
	struct node_s **table;	
};


/* Create a new hashtable. */
struct hashtable_s *createHashTable( int size ) {
	
	struct hashtable_s *newHashTable = NULL;
	int i;
	
	/* Allocate the table itself. */
	newHashTable = malloc( sizeof( struct hashtable_s ) );
	
	/* Allocate pointers to the head nodes. */
	newHashTable->table = malloc( sizeof( struct hashtable_s * ) * size );
	
	/* Setting each index containing a head node to Null*/
	for( i = 0; i < size; i++ ) {
		newHashTable->table[i] = NULL;
	}
	
	newHashTable -> size = size;
	
	return newHashTable;
}

/*Create new node to be entered into hashTable*/
struct node_s *createNode(struct kmer_s val){
	struct node_s *newNode;
	newNode = malloc(sizeof(struct node_s));
	newNode->data=val;
	newNode->next = NULL;
	return newNode;
}

/*Adding a node into a hash table. Utilizes @createNode*/
void addNode(struct hashtable_s *targetTable,int location,struct kmer_s import){
	struct node_s * newEntry = NULL;
	struct node_s * current = NULL;
	struct node_s * previous = NULL;
	newEntry = createNode(import);
	current = targetTable -> table[location];
	
	while( current != NULL){
		previous = current;
		current = current->next;
	}
	if (current == targetTable->table[location]){ //Add at the beginning of a list
		newEntry->next = current;
		targetTable->table[location] = newEntry;
	}else{ //Add at the end of a list
	 previous->next = newEntry;
	}	
}

void removeNode(struct hashtable_s *targetTable,int location,struct kmer_s import){
	struct node_s * toRemove = NULL;
	struct node_s * current = NULL;
	struct node_s * previous = NULL;
	toRemove = createNode(import); //create an instance of what we WANT TO REMOVE
	current = targetTable->table[location];
	while(current!=NULL){
		if((strcmp(current->data.kmerChars,toRemove->data.kmerChars)==0)){
			if(previous == NULL){ //At the start of the list
			 targetTable->table[location] = current -> next;
			return;
			}else{ // Any where else in the list
				 previous->next = current -> next;
			return;
			}
		}
	previous = current;
	current = current -> next;
	}
	printf("No matches found\n");
	return;
}


void printTable(struct hashtable_s *targetTable){
	int i;
	for(i=0;i<targetTable->size;i++){
		struct node_s * current = NULL;
		current = targetTable -> table[i];
		printf("Index i: %i\n", i);
		while(current != NULL){
			printf("%s\n",current->data.kmerChars);
			current = current->next;
		}
	}
}

//Quadratic hashFunction. This will be updated.
int hashCode(char charArray[]){
	int stringLength=strlen(charArray);
	int i;
	int hashNumber = 0;
	for(i=0;i<stringLength;i++){
		char letter = charArray[i];
		int exponentiate;
		exponentiate = pow(letter,i);
		if(i%2==0){
			hashNumber = hashNumber + exponentiate;
		}else{
			hashNumber = hashNumber - exponentiate;
		}
	}
	if(hashNumber < 0){
		hashNumber*=-1;
	}
	return hashNumber;
}

//After creating the hashtable_s struct, use this method with a file to create a healthy hashtable
void populateTable(struct hashtable_s *freshTable, char *sickFileName){
	FILE *fp = fopen(sickFileName,"r");	
	char space[100];								
	while(fgets(space,200,fp)!=NULL){	
		char temp[20];
		int holder = 0;
		memset(&temp[0], 0, sizeof(temp));
		while(space[holder]!=' '){
			temp[holder]=space[holder];
			holder++;
		}
		struct kmer_s newKmer;
		strcpy(newKmer.kmerChars,temp);
		int index = hashCode(newKmer.kmerChars);
		index = index % freshTable->size;
		addNode(freshTable,index,newKmer);	
	}	
}

//Getting rid of all the kmers that are in both the healthy querry and sick querry.
void locateUnique(struct hashtable_s *targetTable, char *healthyFileName){
	FILE *fp = fopen(healthyFileName,"r");	
	char space[100];								
	while(fgets(space,200,fp)!=NULL){	
		char temp[20];
		int holder = 0;
		memset(&temp[0], 0, sizeof(temp));
		while(space[holder]!=' '){
			temp[holder]=space[holder];
			holder++;
		}
		struct kmer_s dummyKmer;
		strcpy(dummyKmer.kmerChars,temp);
		int index = hashCode(dummyKmer.kmerChars);
		index = index % targetTable->size;
		removeNode(targetTable,index,dummyKmer);	
	}
}


void searchThroughFasta(struct hashtable_s *populatedTable,char *fastaFileName,int fileNumber){
	FILE *fp = fopen(fastaFileName,"r");
	char fileOutPutName[20];
	sprintf(fileOutPutName, "file_out%d.txt", fileNumber);
	FILE *fileOutput = fopen(fileOutPutName,"w");
	int fastaCounter=1;
	int DNAid;
	int matchesFound = 0;
	char space[500];
	while(fgets(space,500,fp)!=NULL){
		if(space[0]== '>'){
			DNAid=atoi(space+1);
			continue;
		}
		int i;
		for(i=0;i<populatedTable->size;i++){
			int foundKmer = 0;
			struct node_s * current = NULL;
			current = populatedTable -> table[i];
			while(current != NULL){
				if((strstr(space, current->data.kmerChars)) == NULL){
				current = current->next;
				}else{
					printf("kmer '%s' was found in DNAid %i:\n",current->data.kmerChars,DNAid);
					puts("");
					printf("%s\n",space);
					fprintf(fileOutput,">%i\n",fastaCounter);
					fprintf(fileOutput,"%s",space);
					puts("");
					foundKmer = 1;
					fastaCounter++;
					break;
				}
			}
			if (foundKmer == 1){
				matchesFound++;
				break;
			}
		}
	}
	if(matchesFound == 0){
		puts("No matches were found");
	}
}

int sequenceCounter (FILE *tmpfp){
	FILE *counterPointer = tmpfp;
	int count=0;
	char space[500];
	while(fgets(space,500,counterPointer)!=NULL){
		count++;
	}
	return count;
}

void fileReader(FILE *tmpfp){
	FILE *fileReaderPointer = tmpfp;
	char space[500];
	while(fgets(space,500,fileReaderPointer)!=NULL){
		printf("%s\n",space);
	}
}

FILE *openForSplitFiles(int filecounter) {
	char fileOutPutName[20];
	sprintf(fileOutPutName, "file_part%d.txt", filecounter);
	return fopen(fileOutPutName, "w+");
}

//Creates n number of new files by splitting the fasta files by the entered number of cores used.
void splitFile(char *fastaFileName,int numberOfSplits){
	FILE *fastaFilePtr = fopen(fastaFileName,"r");
	char space[500];
	int numberOfSequences = sequenceCounter(fastaFilePtr);
	rewind(fastaFilePtr);

	int numberOfCores = numberOfSplits;
	int filecounter=1, Sequencecounter=0;
	int split = numberOfSequences/numberOfCores;
	int splitRemainder = numberOfSequences%numberOfCores;
	
	printf("number of sequences: %i\n",numberOfSequences);
	printf("number of cores: %i\n",numberOfCores);
	printf("split numbers: %i %i\n",split,splitRemainder);
	
	FILE *splitFilePtr = openForSplitFiles(filecounter);
	
	while(fgets(space,500,fastaFilePtr)!=NULL){
		if(Sequencecounter != split + splitRemainder){
			fprintf(splitFilePtr,"%s",space);
			Sequencecounter++;
		}else{
			Sequencecounter = 0;
			splitRemainder = 0;
			filecounter++;
			splitFilePtr = openForSplitFiles(filecounter);
			fprintf(splitFilePtr,"%s",space);
		}
	}
}

//These arguements are utilized in each thread
struct threadArgValues{
	char readFileName[20];
	struct hashtable_s *useTable;
	int fileNumber;
};

//The function each thread executes
void* perform_work( void* argument ){
  	struct threadArgValues * instance = argument;
	printf("%s\n",instance->readFileName);
 	 FILE *openForReading = fopen(instance->readFileName,"r");
  	searchThroughFasta(instance->useTable,instance->readFileName,instance->fileNumber);
 	return NULL;
}

void parallelProcessing(int splitCount, struct hashtable_s *thisTable){
	pthread_t threads[splitCount];
  	struct threadArgValues thread_args[splitCount];
  	int result_code;
 	unsigned index;
 	
 	// Creating thread objects to be used in each thread.
 	for(index = 0;index<splitCount;++index){
 		char testReadFileName[15];
 		sprintf(testReadFileName,"file_part%d.txt",index+1);
 		strcpy(thread_args[index].readFileName,testReadFileName);
 		thread_args[index].useTable=thisTable;
 		thread_args[index].fileNumber=index+1;
	}
 	
  	// create all threads one by one
  	for(index = 0;index<splitCount;++index){
    printf("In main: creating thread %d\n", index);
    result_code = pthread_create(threads+index, NULL, perform_work, thread_args+index);
    assert(!result_code);
  	}
  	
  	for(index = 0;index<splitCount;++index){
    result_code = pthread_join(threads[index],NULL);
    assert(!result_code);
  	}
   	exit(EXIT_SUCCESS);
}

main(){

}
