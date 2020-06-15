
# Each cluster is a list of indices x+n*y
# clusters = {}
# For each point, I should have a point to the beginning of a cluster list
# A cluster list: Dictionary
#   cluster_ptr[i] = pointer to cluster list that contains the index i (grid point)
#   cluster list:   header: cluster number, + list of indices
#   list[indx] = cluster_number
#   list[nodes] = list of nodes
#

import numpy as np
import matplotlib.pyplot as plt
import string as s
from sets import Set
import random
import scipy.spatial as kd
from IPython import embed

#import timing

nb_clusters = 0

#t_process = timing.Timer("process")
#t_update = timing.Timer("update")
#t_intersect = timing.Timer("intersect")
#t_chains = timing.Timer("chains")
#t_cluster_init = timing.Timer("cluster init")
#t_tmp = timing.Timer("tmp")
#t_tmp1 = timing.Timer("tmp1")
#t_tmp2 = timing.Timer("tmp2")
#t_kdtree = timing.Timer("kdtree")
#t_kd_nei = timing.Timer("neighbor search")

#----------------------------------------------------------------------
class Set(set):

    def __init__(self, l=None):
        #self.t_set = timing.Timer("Set update")
        #self.t_set.reset()
        if l == None:
            set.__init__(self)
        else:
            set.__init__(self, l)

    #def __init__(self, s):
        #set.__init__(self, s)

    def update(self, s):
        #t_update.tic()
        #self.t_set.tic()
        set.update(self, s)
        #self.t_set.toc()
        #t_update.toc()

    def intersection(self, s):
        #t_intersect.tic()
        r = set.intersection(self, s)
        #print("inside intersection")
        #t_intersect.toc()
        return(r)

    def printUpdateTime(self, msg=""):
        if (msg != ""):
            print("%s" % msg)
        self.t_set.xprint()
        #print "(%s), time: %f sec" % ("Set update", self.t_set.running_time)

#----------------------------
class Cluster:
    """ Efficient kd-tree based cluster algorithm
        Author: Gordon Erlebacher
    """

    which_cluster = {}

    #--------------------------------
    def __init__(self):
        self.max_cluster = 0
        pass

    #--------------------------------
    def __newCluster__(self):
        dictx = {}
        dictx['cluster_nb'] = self.max_cluster
        self.which_cluster[self.max_cluster] = dictx
        self.max_cluster += 1
        dictx['list'] = set([])
        dictx['history'] = set([])  # past incarnations, i.e., cluster addresses
        return(dictx)

    #--------------------------------
    def checkPartioning(self):
        # only use if doubful of results
        for sa in new_list:
            for sb in new_list:
                if sa == sb: continue
                sc = sa.intersection(sb)
                if (len(sc) != 0):
                    print("sa= ", sa)
                    print("sb= ", sb)
                    print("sc= ", sc)
                    print("INCORRECT PARTITIONING")
                    exit()
                #print(sc)
        print("CORRECT PARTIONING")

    #--------------------------------
    def readFile(self, filename):

        fd = open(filename, 'ra')
        array = [[int(x) for x in line.split()] for line in fd]

        x = []
        y = []

        for a in array:
            try:
                x.append(a[0])
                y.append(a[1])
            except:
                break

        arr = np.array( zip(x,y) )
        kdt = kd.cKDTree(arr)
        self.neighbors = kdt.query_ball_point(arr,2)

        return np.array((x,y))

    #------------------------------------
    def __chains__(self):
        set_list = []
        print("__chains__, master_cluster_list len: ", len(self.master_cluster_list))
        for c in self.master_cluster_list:
            cnb = c['cluster_nb']
            h = c['history']
            if len(h) == 0:
                continue
            h.add(cnb)
            set_list.append(h)

        for i in xrange(len(set_list)):
            si = set_list[i]
            for j in xrange(i):
                sj = set_list[j]
                ss = si.intersection(sj)
                if len(ss) > 0:
                    set_list[i].update(sj)
                    set_list[j] = set()

        return set_list

    #----------------------------------
    def pltGrid(self, list, targetShape, grads ):

        pts = np.asarray(targetShape)
        maxV = np.max(grads)
        minV = np.min(grads)
        scale = 1.0/(maxV-minV)
        print "max = %f, min = %f, scale = %f" % (maxV,minV,scale)

        for se in list:
            #if len(se) < 50:
            #    continue
            x = []
            y = []
            cols = []
            #col = colors[co % 7]
            r = random.random()
            g = random.random()
            b = random.random()
            for xy in se:
                #s = grads[tuple(pts[:,xy])]*scale
                s = grads[tuple(pts[:,xy])]
                x.append(pts[0][xy])
                y.append(pts[1][xy])
                cols.append( (r,g,b,s) )

            a = plt.scatter( x, y, color=cols )
            #a = plt.plot(x,y, 'ok' )
            #plt.setp(a, color=col, markersize = 3)

        plt.show()
        return a

#-------------------------------------------------
    def pltGrid1(self, list, pts ):
        # pts: [[x1,x2,...], [y1,y2,...]]
        #colors = ['ob','og','or','oc','om','oy','ok']
        colors = ['r','g','r','c','m','y','k']
        co = 0
        for se in list:
            if len(se) < 50:
                continue
            x = []
            y = []
            col = colors[co % 7]
            #print(se)
            for xy in se:
                x.append(pts[0][xy])
                y.append(pts[1][xy])
            a = plt.plot(x,y, 'ok' )
            plt.setp(a, color=col, markersize = 3)
            co += 1

        plt.show()
    #-----------------------------------------------
    def expand(self, pt_nb):
        count = 0
        seeds = []
        nei = self.neighbors[pt_nb]
        #print("nb nei: ", len(nei))
        source_cluster_ptr = self.cluster_ptr[pt_nb]
        # neighbors contain the center (at an arbitrary location)
        for i in nei:

            if self.cluster_ptr[i] != None:  #visited
                continue

            seeds.append(i)
            count += 1
            #if (count % 10): print count
            self.cluster_ptr[i] = source_cluster_ptr

        while (len(seeds) > 0):
            #print "len(seeds)= ", len(seeds)
            seed = seeds.pop()
            #print "seed: ", seed

            #if self.cluster_ptr[seed] != None:
                #continue

            seed_neighbors = self.neighbors[seed]
            #print "len(seed_neighbors): ", len(seed_neighbors)

            for s in seed_neighbors:
                # seed has already has a cluster number
                if self.cluster_ptr[s] == None:
                    seeds.append(s)
                    count += 1
                    #if (count % 10): print count
                    self.cluster_ptr[s] = source_cluster_ptr;

        #print("nb seeds= ", count)
        #exit()
    """
       seeds = neighbors of pt_nb (excluding pt_nb)
       while (len(seeds) > 0):
         seed = seeds.pop()
         seeds += neighbors(seed)  # for seeds without classification
       }
    """
    #-----------------------------------------------
    def __processPoint__(self, pt_nb):
        global nb_clusters

        # if pt_nb has been classfied, exit
        # assign new cluster to pt

        if self.cluster_ptr[pt_nb] == None:
            self.cluster_ptr[pt_nb] = self.__newCluster__()
            nb_clusters += 1
            clust_nb = nb_clusters-1
            self.cluster_ptr[pt_nb]['cluster_nb'] = clust_nb
            self.expand(pt_nb)
            #print "new cluster at point: ", pt_nb

        return

    #-----------------------------------------------
    def clusterGrid( self, pts, radius ):

        #print("\nClusterGrid")

        #t_cluster_init.tic()
        try:
            self.neighbors
        except:
            arr = np.array(pts.transpose())
            #t_kdtree.tic()
            kdt = kd.cKDTree(arr)
            #t_kdtree.toc(out=1)
            #t_kd_nei.tic()
            self.neighbors = kdt.query_ball_point(arr,radius)
            #t_kd_nei.toc(out=1)

        self.max_cluster = 0

        clusters = []
        self.master_cluster_list = []

        nb_pts = pts.shape[1]

        clusters = []
        self.cluster_ptr = [None]*nb_pts
        self.master_cluster_list = []
        #t_cluster_init.toc(out=1)

        #t_process.tic()
        for pc in xrange(nb_pts):
            self.__processPoint__(pc)  # pts is list of (x,y)
        #t_process.toc(out=1)

        #t_tmp.xprint()
        #t_tmp1.xprint()
        #t_tmp2.xprint()



        ## create cluster lists
        #print("len master_cluster_list: ", len(self.master_cluster_list))
        #print("nb clusters: ", nb_clusters)
        new_list = [0] * nb_clusters
        for pc in xrange(nb_pts):
            cnb = self.cluster_ptr[pc]['cluster_nb']
            try:
                new_list[cnb].append(pc)
            except:
                new_list[cnb] = []
                new_list[cnb].append(pc)

        return new_list


        #t_chains.tic()
        chain = self.__chains__()
        all_chain_indices = set()
        for c in chain:
            all_chain_indices.update(c)

        all_cluster_indices = set()
        for c in self.master_cluster_list:
            all_cluster_indices.update([c['cluster_nb']])

        remaining_indices = all_cluster_indices.difference(all_chain_indices)

        # link the lists of each chain
        old_list = self.master_cluster_list
        new_list = []

        # Construct master list from chains
        for ch in chain:
            el = set()
            for c in ch:
                #print(self.master_cluster_list[c]['list'])
                el.update(self.master_cluster_list[c]['list'])
            if (len(el) == 0):
                continue
            new_list.append(el)

        # Add the remaining elements (those not in all_indices)

        for c in remaining_indices:
            new_list.append(set(self.master_cluster_list[c]['list']))

        #t_chains.toc(out=1)

        return new_list

#----------------------------------------------------------------------
if __name__ == "__main__":

    cluster = Cluster()
    pts = cluster.readFile("data/ocean.dat")

    cluster_list = cluster.clusterGrid(pts)

    print("press enter")
    raw_input()
    cluster.pltGrid(cluster_list, pts)

#----------------------------------------------------------------------

