#!/usr/bin/env python3
### 2023.08.07
### a script for building a neural network for identifying soft-clipped read pairs
### 2023.09.03
### a mobilenetV3 model is added
### 2023.09.25
### train with simulated data


import numpy as np
import torch
from torch.utils.data import Dataset
from torch.utils.data import DataLoader
import torch.nn as nn
import torch.nn.functional as F


### construct very simple simulation dataset
# import csv
# testfile=open("/home/zhangyf/rd/TE_insertion/SomTD/testGsSim/newSimChr1/testSup/testDeep/test.resource.txt","w")
# tsv_output = csv.writer(testfile, delimiter='\t')
# for i in range(100):
#     for j in range(10):
#         for k in range(36):
#             _=tsv_output.writerow(list(np.random.normal(j,0.1,size=150)))
# testfile.close()

### construct real simulation dataset, only once
from SomTD_ReadGTP import *
import csv
import pysam
trainDataFile=pysam.AlignmentFile("/home/zhangyf/rd/TE_insertion/SomTD/CNNdata/trainData/trainData.combine.clip.pair.sam","r")
trainDataTsvfile=open("/home/zhangyf/rd/TE_insertion/SomTD/CNNdata/trainData/trainData.tsv","w")
tsv_output = csv.writer(trainDataTsvfile, delimiter='\t')
trainDataList=[]
NR = 1
interList=[]
for i in trainDataFile:
    interList.append(i)
    if NR % 12 == 0:
        trainDataList.append(interList)
        interList = []
    NR = NR + 1

trainDataFile.close()
TElenDict = {'ALU':281,'LINE1':6019,'SVA':1316,'HERVK':8759}

for i in trainDataList:
    testGTP=ReadGTP(i,150,TElenDict)
    testGTP.aligns2array()
    _=tsv_output.writerows(testGTP.alignArray)

trainDataTsvfile.close()

# define dataset

class AlignmentDataset(Dataset):
    """
    a dataset class of simulated TE insertion alignments for the model
    """
    def __init__(self,dataFileLoc,labelFileLoc):
        """
        method for constructing alignment array list with a shape:  and alignment array label list
        """
        super(AlignmentDataset, self).__init__()
        dataList=[]
        labelList=[]
        # construct array list
        combineFile = open(dataFileLoc,"r")
        NR = 1
        valueList = [] # variable for saving a single channel
        channelList = [] # variable for saving 18 channels
        for i in combineFile:
            valueList.append(list(map(float,i.strip("\n").split("\t"))))
            if NR % 2 == 0:
                channelList.append(valueList)
                valueList = []
            if NR % 36 == 0:
                dataList.append(torch.tensor(channelList,dtype=torch.float))
                channelList = []
            NR = NR + 1
        # construct label list
        labelFile = open(labelFileLoc,"r")
        for i in labelFile:
            if len(i.strip("\n").split("\t")) ==1:
                labelList.append(torch.tensor(int(i.strip("\n")),dtype=torch.long)) # torch.long type for label
                # labelList.append(torch.cuda.FloatTensor(int(i.strip("\n"))))
            else:
                labelList.append(torch.tensor(int(i.strip("\n").split("\t")[1]),dtype=torch.long)) # torch.long type for label
                # labelList.append(torch.cuda.Tensor(int(i.strip("\n").split("\t")[1])))
        self.dataList = dataList
        self.labelList = labelList
    def __len__(self):
        return len(self.dataList)  # number of samples in the dataset
    def __getitem__(self, index):
        return self.dataList[index], self.labelList[index]



# construct data and dataloader
# trainingData = AlignmentDatasetFC("/home/zhangyf/rd/TE_insertion/SomTD/testGsSim/newSimChr1/testSup/testDeep/test.resource.txt","/home/zhangyf/rd/TE_insertion/SomTD/testGsSim/newSimChr1/testSup/testDeep/test.label.txt")
# trainingData = AlignmentDataset("/home/zhangyf/rd/TE_insertion/SomTD/testGsSim/newSimChr1/testSup/testDeep/test.resource.txt","/home/zhangyf/rd/TE_insertion/SomTD/testGsSim/newSimChr1/testSup/testDeep/test.label.txt")
trainingData = AlignmentDataset(r"D:\OneDrive - mails.ucas.ac.cn\trainData\trainData.tsv",r"D:\OneDrive - mails.ucas.ac.cn\trainData\trainData.index")

# testData = AlignmentDatasetFC()
trainDataloader = DataLoader(trainingData, batch_size=64, shuffle=True)
# testDataloader = DataLoader(trainingData, batch_size=64, shuffle=False)

# define models
class FCNet(nn.Module):
    """
    a two-layer fully connected model for detecting nonreference TE insertion read pairs
    """
    def __init__(self):
        super(FCNet, self).__init__()
        self.fc1 = nn.Linear(5400, 256)
        self.fc2 = nn.Linear(256,128)
        self.fc3 = nn.Linear(128, 10)
        self.dropout1 = nn.Dropout(0.25)
    def forward(self, x):
        x = x.flatten(1) # Flatten x with start_dim=1
        x = self.fc1(x) # hidden layer
        x = F.relu(x)
        x = self.fc2(x)
        x = F.relu(x)
        x = self.dropout1(x)
        x = self.fc3(x) # output layer
        output = F.log_softmax(x, dim=1) # make sure every row instead of column is normalized
        return output

modelFCNet = FCNet()

class LeNet(nn.Module):
    """
    a LENET model for detecting nonreference TE insertion read pairs
    """
    def __init__(self):
        super(LeNet, self).__init__()
        self.conv1 = nn.Conv2d(in_channels = 18, out_channels = 36, kernel_size = (1,3), stride = (1,1), padding = (0,0))
        self.maxpool1 = nn.MaxPool2d(kernel_size = (1,4), stride = (1,2))
        self.conv2 = nn.Conv2d(in_channels = 36, out_channels = 72, kernel_size = (1,3), stride = (1,2), padding = (0,0))
        self.maxpool2 = nn.MaxPool2d(kernel_size = (1,2), stride = (1,2))
        self.fc1 = nn.Linear(72 * 2 * 18, 120)
        self.fc2 = nn.Linear(120, 84)
        self.fc3 = nn.Linear(84, 10)
        self.dropout1 = nn.Dropout(0.25)
        self.dropout2 = nn.Dropout(0.25)
    def forward(self, x): # N * C * H * W initiation: 64 * 18 * 2 * 150 
        x = self.conv1(x) # 64 * 36 * 2 * 
        x = F.relu(x)
        x = self.maxpool1(x) # 64 * 36 * 2 * 36
        x = self.conv2(x) # 64 * 72 * 2 * 17
        x = F.relu(x)
        x = self.maxpool2(x) # 64 * 72 * 2 * 8
        x = x.reshape(x.size()[0], -1) # 64 * 1152
        x = self.fc1(x)
        x = F.relu(x)
        x = self.dropout1(x)
        x = self.fc2(x)
        x = F.relu(x)
        x = self.dropout2(x)
        x = self.fc3(x)
        output = F.log_softmax(x, dim=1) # make sure every row instead of column is normalized
        return output

modelLeNet = LeNet()


class SEModule(nn.Module):
    """
    Squeeze-and-Excite module
    """
    def __init__(self, in_size, reduction=4):
        super(SEModule, self).__init__()
        squeeze_size =  max(in_size // reduction, 8)
        self.se = nn.Sequential(
            nn.AdaptiveAvgPool2d(1),
            nn.Conv2d(in_size, squeeze_size, kernel_size=1, bias=False),
            nn.BatchNorm2d(squeeze_size),
            nn.ReLU(inplace=True),
            nn.Conv2d(squeeze_size, in_size, kernel_size=1, bias=False),
            nn.Hardsigmoid()
        )
    def forward(self, x):
        return x * self.se(x)


class Block(nn.Module):
    """
    Mobilenet V3 block containing a pointwise expansion layer, a depthsise layer, a optional SE
    """
    def __init__(self, kernel_size, in_size, expand_size, out_size, act, se, stride):
        super(Block, self).__init__()
        self.stride = stride
        self.in_size = in_size
        self.expand_size = expand_size
        self.out_size = out_size
        #
        self.conv1 = nn.Conv2d(in_size, expand_size, kernel_size=1, bias=False)
        self.bn1 = nn.BatchNorm2d(expand_size)
        self.act1 = act(inplace=True)
        #
        self.conv2 = nn.Conv2d(expand_size, expand_size, kernel_size=(1,kernel_size), stride=(1,stride), padding=(0,kernel_size//2), groups=expand_size, bias=False)
        self.bn2 = nn.BatchNorm2d(expand_size)
        self.act2 = act(inplace=True)
        self.se = SEModule(expand_size) if se else nn.Identity()
        #
        self.conv3 = nn.Conv2d(expand_size, out_size, kernel_size=1, bias=False)
        self.bn3 = nn.BatchNorm2d(out_size)
        self.act3 = act(inplace=True)
    def forward(self, x):
        skip = x
        #
        out = self.act1(self.bn1(self.conv1(x)))
        out = self.act2(self.bn2(self.conv2(out)))
        out = self.se(out)
        out = self.bn3(self.conv3(out))
        #
        if self.stride == 1 and self.in_size == self.out_size:
            return out + skip
        else:
            return out


class MobileNetV3_Small(nn.Module):
    """
    a small version MobileNetV3 for detecting nonreference TE insertion read pairs
    """
    def __init__(self, num_classes=10, act=nn.Hardswish):
        super(MobileNetV3_Small, self).__init__()
        self.conv1 = nn.Conv2d(18, 18, kernel_size=(1,3), stride=(1,2), padding=(0,1), bias=False)
        self.bn1 = nn.BatchNorm2d(18)
        self.hs1 = act(inplace=True)
        #
        self.bneck = nn.Sequential(
            Block(3, 18, 18, 18, nn.ReLU, True, 2),
            Block(3, 18, 72, 24, nn.ReLU, False, 2),
            Block(3, 24, 88, 24, nn.ReLU, False, 1),
            Block(5, 24, 96, 40, act, True, 2),
            Block(5, 40, 240, 40, act, True, 1),
            Block(5, 40, 240, 40, act, True, 1),
            Block(5, 40, 120, 48, act, True, 1),
            Block(5, 48, 144, 48, act, True, 1),
            Block(5, 48, 288, 96, act, True, 2),
            Block(5, 96, 576, 96, act, True, 1),
            Block(5, 96, 576, 96, act, True, 1),
        )
        #
        self.conv2 = nn.Conv2d(96, 576, kernel_size=1, stride=1, padding=0, bias=False)
        self.bn2 = nn.BatchNorm2d(576)
        self.hs2 = act(inplace=True)
        self.gap = nn.AdaptiveAvgPool2d(1)
        #
        self.linear3 = nn.Linear(576, 1280)
        self.bn3 = nn.BatchNorm1d(1280)
        self.hs3 = act(inplace=True)
        self.drop = nn.Dropout(0.2,inplace=True)
        self.linear4 = nn.Linear(1280, num_classes)
        self.init_params()
    def init_params(self):
        for m in self.modules():
            if isinstance(m, nn.Conv2d):
                nn.init.kaiming_normal_(m.weight, mode='fan_out')
                if m.bias is not None:
                    nn.init.zeros_(m.bias)
            elif isinstance(m, (nn.BatchNorm2d, nn.GroupNorm)):
                nn.init.ones_(m.weight)
                nn.init.zeros_(m.bias)
            elif isinstance(m, nn.Linear):
                nn.init.normal_(m.weight, 0, 0.01)
                nn.init.zeros_(m.bias)
    def forward(self, x):
        out = self.hs1(self.bn1(self.conv1(x)))
        out = self.bneck(out)
        #
        out = self.hs2(self.bn2(self.conv2(out)))
        out = self.gap(out).flatten(1)
        out = self.drop(self.hs3(self.bn3(self.linear3(out))))
        #
        return self.linear4(out)

modelMBN = MobileNetV3_Small(num_classes=10)



# define the loss function and optimizer
# cross entropy and adam are used
# from torch.optim import RMSprop
from torch.optim import Adam
from torch.optim import SGD
lossFunc = nn.CrossEntropyLoss()
optimizer = SGD(modelMBN.parameters(),lr=0.001)
optimizer = Adam(modelMBN.parameters(), lr=0.0001) # performs better on simple data
# optimizer = RMSprop(modelMBN.parameters()) # default parameters performs better than momentum = 0.9, though much more worse than Adam


# define train and test function
def trainModel(model,epochNum):
    """
    function for training the neural network
    """
    device = torch.device("cuda:0" if torch.cuda.is_available() else "cpu") # model check
    print("The model will be trained on", device, "device")
    model.to(device)
    epochLossList=[]
    epochAccList=[]
    for epoch in range(epochNum):  # loop over the dataset multiple times
        epochLoss = 0.0
        epochAcc = 0.0
        endMark = 0
        for i,data in enumerate(trainDataloader): # i is batch ID
            input,label=data
            input = input.to(device)
            label = label.to(device)
            optimizer.zero_grad()
            output=model(input)
            loss=lossFunc(output,label) # generate loss
            epochAcc += (output.argmax(1) == label).type(torch.float).sum().item() # generate accuracy
            loss.backward()
            optimizer.step()
            epochLoss += loss.item()
            if i % 3134 == 0 and i != 0:
                print("epoch " + str(epoch + 1) + " loss: " + str(round(epochLoss/200600,4)) + "; accuracy: " + str(round(epochAcc/200600,4)) + "       ",end="\r")
                epochLossList.append(round(epochLoss/200600,4))
                epochAccList.append(round(epochAcc/200600,4))
                if epochAcc/200600 >= 0.99:
                    endMark = 1
                epochAcc = 0
                epochLoss = 0
        if endMark == 1:
            print("\n" + "epoch " + str(epoch + 1) + " is enough")
            # torch.save(model.state_dict(), "/home/zhangyf/rd/TE_insertion/SomTD/testGsSim/newSimChr1/testSup/testDeep/modelMBN.pth") # save model
            torch.save(model.state_dict(), r"D:\OneDrive - mails.ucas.ac.cn\trainData\modelMBN.pth")
            break
    return epochLossList,epochAccList

train1,train2=trainModel(modelMBN,100)
num_params = sum(p.numel() for p in modelMBN.parameters())


# load model
modelMBN = MobileNetV3_Small()
modelMBN.load_state_dict(torch.load("/home/zhangyf/rd/TE_insertion/SomTD/testGsSim/newSimChr1/testSup/modelMBN.pth"))
modelMBN.eval()
modelMBN(torch.tensor(np.random.normal(1,0.1,size=(1,18,2,150)),dtype=torch.float))


### 2023.10.06 training results
# torch.save(modelMBN.state_dict(), r"D:\OneDrive - mails.ucas.ac.cn\trainData\modelMBN.pth")
train1=[0.0206, 0.0167, 0.0143, 0.0089, 0.0068, 0.0059, 0.0053, 0.005, 0.0047, 0.0045, 0.0043, 0.0042, 0.0041, 0.004, 0.0039, 0.0038, 0.0037, 0.0036, 0.0035, 0.0035, 0.0034, 0.0033, 0.0032, 0.0032, 0.0031, 0.0031, 0.003, 0.003, 0.0029, 0.0028, 0.0028, 0.0027, 0.0027, 0.0026, 0.0026, 0.0025, 0.0025, 0.0024, 0.0024, 0.0023, 0.0023, 0.0023, 0.0022, 0.0022, 0.0021, 0.0021, 0.0021, 0.002, 0.002, 0.0019, 0.0019, 0.0019, 0.0018, 0.0018, 0.0018, 0.0018, 0.0017, 0.0017, 0.0017, 0.0016, 0.0016, 0.0016, 0.0016, 0.0015, 0.0015, 0.0015, 0.0015, 0.0014, 0.0014, 0.0014, 0.0014, 0.0013, 0.0013, 0.0013, 0.0013, 0.0013, 0.0013, 0.0013, 0.0012, 0.0012, 0.0012, 0.0012, 0.0012, 0.0012, 0.0011, 0.0011, 0.0011, 0.0011, 0.0011, 0.0011, 0.0011, 0.0011, 0.0011, 0.001, 0.0011, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.0009, 0.001, 0.0009, 0.0009, 0.0009, 0.0009, 0.0009, 0.0009, 0.0009, 0.0009, 0.0009, 0.0009, 0.0009, 0.0008, 0.0009, 0.0008, 0.0008, 0.0008, 0.0008, 0.0008, 0.0008, 0.0008, 0.0008, 0.0008, 0.0008, 0.0008, 0.0008, 0.0008, 0.0008, 0.0008, 0.0007, 0.0008, 0.0007, 0.0007, 0.0008, 0.0007, 0.0007, 0.0007, 0.0007, 0.0007, 0.0007, 0.0007, 0.0007, 0.0007, 0.0007, 0.0007, 0.0007, 0.0006, 0.0007, 0.0006, 0.0007, 0.0007, 0.0006, 0.0006, 0.0006, 0.0006, 0.0006, 0.0006, 0.0006, 0.0006, 0.0006, 0.0006, 0.0006, 0.0006, 0.0006, 0.0006, 0.0006, 0.0006, 0.0006, 0.0006, 0.0006, 0.0006, 0.0006, 0.0005, 0.0006, 0.0005, 0.0005, 0.0005, 0.0005, 0.0005, 0.0005, 0.0005, 0.0005, 0.0005, 0.0005, 0.0005, 0.0005, 0.0005, 0.0005, 0.0005, 0.0005, 0.0005, 0.0005, 0.0005, 0.0005, 0.0005]

train2=[0.4267, 0.4849, 0.5649, 0.7671, 0.8281, 0.8547, 0.8693, 0.878, 0.8856, 0.89, 0.8945, 0.8973, 0.9005, 0.9033, 0.9055, 0.9079, 0.9096, 0.9114, 0.913, 0.9147, 0.9169, 0.9177, 0.9196, 0.9209, 0.9234, 0.9242, 0.9255, 0.9259, 0.9283, 0.9296, 0.9305, 0.9313, 0.9329, 0.9332, 0.935, 0.9361, 0.937, 0.9379, 0.9394, 0.9406, 0.9404, 0.9421, 0.9436, 0.9442, 0.9457, 0.9462, 0.9476, 0.9491, 0.9495, 0.95, 0.9511, 0.9526, 0.9524, 0.9535, 0.9535, 0.9543, 0.9559, 0.9569, 0.9565, 0.9578, 0.9579, 0.9586, 0.9597, 0.9604, 0.9614, 0.9614, 0.962, 0.9625, 0.9635, 0.9638, 0.9637, 0.9647, 0.9652, 0.9654, 0.9661, 0.9666, 0.9668, 0.9675, 0.9674, 0.968, 0.9688, 0.9682, 0.9692, 0.9694, 0.9698, 0.9697, 0.971, 0.9703, 0.9711, 0.9709, 0.9715, 0.9714, 0.972, 0.9722, 0.9721, 0.9725, 0.9727, 0.973, 0.9736, 0.9741, 0.974, 0.9747, 0.9748, 0.9749, 0.9747, 0.9746, 0.9751, 0.9753, 0.9758, 0.9758, 0.9757, 0.9764, 0.9766, 0.9765, 0.9762, 0.9771, 0.9772, 0.9778, 0.9775, 0.9781, 0.9781, 0.9777, 0.9781, 0.9782, 0.9786, 0.9788, 0.9792, 0.9791, 0.9789, 0.9796, 0.9795, 0.9802, 0.9793, 0.9801, 0.9806, 0.9801, 0.9802, 0.9805, 0.9805, 0.9812, 0.9803, 0.9813, 0.9811, 0.981, 0.9819, 0.9815, 0.9818, 0.9818, 0.9825, 0.9821, 0.9824, 0.9832, 0.9822, 0.9833, 0.9832, 0.9833, 0.9832, 0.9831, 0.9839, 0.9833, 0.9843, 0.9833, 0.9841, 0.9845, 0.9843, 0.9842, 0.9847, 0.9843, 0.9849, 0.9845, 0.9853, 0.9849, 0.9853, 0.9851, 0.9856, 0.9854, 0.9857, 0.9863, 0.9852, 0.9861, 0.9862, 0.9865, 0.9863, 0.9867, 0.9865, 0.9867, 0.987, 0.9869, 0.9871, 0.9873, 0.9873, 0.9873, 0.9874, 0.9875, 0.9874, 0.988, 0.9878, 0.9882, 0.9882, 0.9881]


def testModel():
    """
    function of testing the neural network
    """