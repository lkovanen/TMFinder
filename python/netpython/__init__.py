""" A module containing networking framework and analysis functions

>>> from netpython import *
>>> def basicTests(net):
...  print len(net.edges)
...  print net[1].deg()
...  print net[1].strength()
...  print percolator.getComponents(net).getSizeDist()
...  print sum(net.weights)
...  #print netext.mst(net).edges
...  #print netext.deg(net)
...  #print netext.clustering(net)
>>> def addNodes1(net):
...  net[1,2]=2.0
...  net[2,1]=1.0
...  net[2,3]=2.0
...  net[3,1]=3.0
...  net[3,4]=4.0
...  net[5,6]=5.0
>>> def addNodes2(net):
...  net['foo','bar']=10.0
...  net['bar','foo']=20.0
...  net['cat','dog']=30.0
>>> net1=pynet.SymmNet()
>>> addNodes1(net1)
>>> basicTests(net1)
5
2
4.0
{2: 1, 4: 1}
15.0
>>> net2=pynet.Net()
>>> addNodes1(net2)
>>> basicTests(net2)
6
1
2.0
{2: 1, 4: 1}
17.0
>>> net3=pynet.SymmFullNet(10)
>>> addNodes1(net3)
>>> basicTests(net3)
5
2
4.0
{2: 1, 4: 1}
15.0
>>> net4=pynet.FullNet(10)
>>> addNodes1(net4)
>>> basicTests(net4)
6
1
2.0
{2: 1, 4: 1}
17.0
>>> net5=pynet.SymmNet()
>>> addNodes1(net5)
>>> addNodes2(net5)
>>> basicTests(net5)
7
2
4.0
{2: 3, 4: 1}
65.0
"""
__all__=["netio","pynet","netext","percolator","eden","visuals","netanalysis","models","transforms","communities","dialogues"]
