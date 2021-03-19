import numpy as np
from collections import OrderedDict

class Lab2(object):
    
    def smith_waterman_alignment(self,s1,s2,penalties) :
        '''
        Input - two sequences and a dictionary with penalities for match, mismatch and gap
        Output - an integer value which is the maximum smith waterman alignment score
        '''
        #start code here
        len_x=len(s1)
        len_y=len(s2)
        dp=np.zeros((len_y+1,len_x+1))
        for i in range(len_x+1):
            dp[0][i]=0
        for j in range(len_y+1):
            dp[j][0]=0
        for i in range(1,len_x+1):
            for j in range(1,len_y+1):
                if s1[i-1]==s2[j-1]:
                    dp1=dp[j-1][i-1]+penalties['match']
                else:
                    dp1=dp[j-1][i-1]+penalties['mismatch']
                dp[j][i]=max(dp[j-1][i]+penalties['gap'],dp[j][i-1]+penalties['gap'],dp1,0)
       
        return int(dp.max())
        #end code here

    def print_smith_waterman_alignment(self,s1,s2,penalties) :
        '''
        Input - two sequences and a dictionary with penalities for match, mismatch and gap
        Output - a tuple with two strings showing the two sequences with '-' representing the gaps
        '''
        #start code here
        max_score=0
        len_x=len(s1)
        len_y=len(s2)
        dp=np.zeros((len_x+1,len_y+1))
        for i in range(len_x+1):
            dp[i][0]=0
        for j in range(len_y+1):
            dp[0][j]=0
        for i in range(len_x):
            for j in range(len_y):
                if s1[i]==s2[j]:
                    dp1=dp[i][j]+penalties['match']
                else:
                    dp1=dp[i][j]+penalties['mismatch']
                temp=max(dp[i][j+1]+penalties['gap'],dp[i+1][j]+penalties['gap'],dp1,0)
                if temp>max_score:
                    max_score=temp
                    max_i = i+1
                    max_j = j+1
                dp[i+1][j+1] = temp
       
        i = max_i
        j = max_j
        output1 = ""
        output2 = ""
        while dp[i][j]:
            if s1[i-1] == s2[j-1]:
                output1 = s1[i-1]+output1
                output2 = s2[j-1]+output2
                i -= 1
                j -= 1
            elif dp[i][j-1]+penalties['gap'] == dp[i][j]:
                output1 = "-"+output1
                output2 = s2[j-1]+output2
                j -= 1
            elif dp[i-1][j]+penalties['gap'] == dp[i][j]:
                output1 = s1[i-1]+output1
                output2 = "-"+output2
                i -= 1
            else:
                output1 = s1[i-1]+output1
                output2 = s2[j-1]+output2
                i -= 1
                j -= 1
                
        return (output1, output2)
                
        
        #end code here

    def find_exact_matches(self,list_of_reads,genome):
       
        '''
        Input - list of reads of the same length and a genome fasta file (converted into a single string)
        Output - a list with the same length as list_of_reads, where the ith element is a list of all locations (starting positions) in the genome where the ith read appears. The starting positions should be specified using the "chr2:120000" format
        '''
        raw= genome.split('>')[1:]
        data=[]
        gen_dict={}
        result=[]
        number=['0','1','2','3','4','5','6','7','8','9']
        for rd in raw:
            rd=rd.replace('\n','')
            for e in number:
                rd=rd.replace(e,'')
            data.append(rd.replace('chr',''))
        #print(data[0:2])
        length=len(list_of_reads[0])
        for i in range(len(data)):
            if len(data[i])>=length:
                for j in range(len(data[i])-length+1):
                    if data[i][j:j+length] not in gen_dict:
                        label='chr'+str(i+1)+':'+str(j+1)
                        gen_dict[data[i][j:j+length]]=[label]
                    else:
                        label='chr'+str(i+1)+':'+str(j+1)
                        gen_dict[data[i][j:j+length]].append(label)
        #print(gen_dict)
        for e in list_of_reads:
            if e in gen_dict:
                result.append(gen_dict[e])
            else:
                result.append([])
                
        return result
        
        #start code here
        #end code here
       
    
    def find_approximate_matches(self,list_of_reads,genome):
        '''
        Input - list of reads of the same length and a genome fasta file (converted into a single string)
        Output -  a list with the same length as list_of_reads, where the ith element is a list of all locations (starting positions) in the genome which have the highest smith waterman alignment score with ith read in list_of_reads
        '''
        raw= genome.split('>')[1:]
        data=[]
        gen_dict={}
        result=[]
        max_score=0
        max_list=[]
        number=['0','1','2','3','4','5','6','7','8','9']
        penalties={'match':1,'mismatch':-1,'gap':-1}
        for rd in raw:
            rd=rd.replace('\n','')
            for e in number:
                rd=rd.replace(e,'')
            data.append(rd.replace('chr',''))
        
        #print(data)
        length=int(len(list_of_reads[0])/4)
        #print(length)
        for i in range(len(data)):
            if len(data[i])>=length:
                for j in range(len(data[i])-length+1):
                    if data[i][j:j+length] not in gen_dict:
                        label='chr'+str(i+1)+':'+str(j+1)
                        gen_dict[data[i][j:j+length]]=[label]
                    else:
                        label='chr'+str(i+1)+':'+str(j+1)
                        gen_dict[data[i][j:j+length]].append(label)
        map=[[0]*(len(list_of_reads[0])-length+1) for _ in range(len(list_of_reads))]
        for i in range(len(list_of_reads)):
            for j in range(len(list_of_reads[0])-length+1):
                segment=list_of_reads[i][j:int(j+len(list_of_reads[0])/4)]
                if segment in gen_dict:
                    map[i][j]=gen_dict[segment]
        #score=[[0]*len(list_of_reads[0])-length+1 for _ in range(len(list_of_reads))]
        #print(map[1][2])
        #print(map)
        for i in range(len(list_of_reads)):
            max_list=[]
            for j in range(len(map[i])):
                if map[i][j]!=0:
                    for k in range(len(map[i][j])):
                        msg=map[i][j][k]
                        channel=int(msg.split(':')[0].split('chr')[1])
                        idx=int(msg.split(':')[-1])
                        x=list_of_reads[i]
                        if idx-j-1+len(x)<=len(data[channel-1]):
                            y=data[channel-1][idx-1-j:idx-j-1+len(x)]
                       # print(len(x)==len(y))
                            score=self.smith_waterman_alignment(x,y,penalties)
                        
                            if len(max_list)==0:
                                max_score=score
                                if ('chr'+str(channel)+':'+str(idx-j)) not in max_list:
                                    max_list.append('chr'+str(channel)+':'+str(idx-j))
                            else:
                                if len(max_list)>0 and score==max_score:
                                    if ('chr'+str(channel)+':'+str(idx-j)) not in max_list:
                                        max_list.append('chr'+str(channel)+':'+str(idx-j))
                                elif score>max_score:
                                    max_score=score
                                    max_list=[]
                                    if ('chr'+str(channel)+':'+str(idx-j)) not in max_list:
                                        max_list.append('chr'+str(channel)+':'+str(idx-j))
                                
                        
                        
#                         if len(max_list)>0 and score==max_score:
#                             max_list.append('chr'+str(channel)+':'+str(idx-j))
#                         elif len(max_list)==0:
#                             max_score=score
#                             max_list.append('chr'+str(channel)+':'+str(idx-j))
#                         else :
#                             if score>max_score:
#                                 max_score=score
#                                 max_list=[]
#                                 max_list.append('chr'+str(channel)+':'+str(idx-j))
                #print(max_list)
            result.append(max_list)   
        return result

            
        
        #start code here
        #end code here
        