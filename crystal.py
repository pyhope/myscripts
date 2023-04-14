import os
import numpy as np
import crystal
import re
import copy


class crystal:
    # 晶胞类,属性和方法见本文件中的help函数

    def __init__(self):
        self.name = "example_MgO"
        # 定义晶胞名称
        self.cell_parameter = {"a": [4.2109999657, 0.0000000000, 0.0000000000], "b": [
            0.0000000000, 4.2109999657, 0.0000000000], "c": [0.0000000000, 0.0000000000, 4.2109999657]}
        # 定义晶胞参数
        self.atomtype = {"Mg": 4, "O": 4}
        # 定义晶胞中原子总类和每种的数目
        self.atomtypelist = ["Mg", "O"]
        # 定义晶胞中原子的种类列表
        self.position = {("Mg", 0): [0.000000000, 0.000000000, 0.000000000], ("Mg", 1): [0.000000000, 0.500000000, 0.500000000], ("Mg", 2): [0.500000000, 0.000000000, 0.500000000], ("Mg", 3): [0.500000000, 0.500000000, 0.000000000], ("O", 0): [
            0.500000000, 0.500000000, 0.500000000], ("O", 1): [0.500000000, 0.000000000, 0.000000000], ("O", 2): [0.000000000, 0.500000000, 0.000000000], ("O", 3): [0.000000000, 0.000000000, 0.500000000]}
        # 定义晶胞中各原子位置字典，key为原子名+序号，序号以0开始。
        # 以MgO为例定义晶胞

    def __str__(self):
        return "晶胞:%s|总原子数:%d|是否自洽:%d" % (self.name, len(self.position), self.check())

    def atomdistance(self, left, right):
        if (left in self.position) and (right in self.position):
            left_array = self.posvector(left)
            cell_array = self.CellParameter()
            for x in [-1.0, 0.0, 1.0]:
                for y in [-1.0, 0.0, 1.0]:
                    for z in [-1.0, 0.0, 1.0]:
                        right_array = self.posvector(right)
                        right_array = right_array+np.array([x, y, z])
                        cell_vector = right_array-left_array
                        real_vector = np.dot(cell_vector, cell_array)
                        dis = np.linalg.norm(real_vector)
                        try:
                            dislist.append(dis)
                        except:
                            dislist = [dis]
            return min(dislist)
        else:
            return -1

    def CellParameter(self):
        return np.array([self.cell_parameter["a"], self.cell_parameter["b"], self.cell_parameter["c"]])

    def check(self):
        # 对晶胞进行自洽检查，包括：
        # 1.原子总数是否与位置字典中key数匹配，否则返回0
        atomtypes = 0
        for atom in self.atomtype.keys():
            atomtypes += self.atomtype[atom]
        if atomtypes != len(self.position):
            return 0
        else:
            return 1

    def clear(self):
        self.atomtype.clear()
        self.atomtypelist.clear()
        self.position.clear()
        self.cell_parameter = {"a": [0.0, 0.0, 0.0], "b": [
            0.0, 0.0, 0.0], "c": [0.0, 0.0, 0.0]}
        self.name = "empty_cell"

    def copy(self):
        out = copy.deepcopy(self)
        out.position = self.position.copy()
        out.atomtype = self.atomtype.copy()
        out.cell_parameter = self.cell_parameter.copy()
        out.atomtypelist = self.atomtypelist.copy()
        return out

    def expend(self, xyz):
        if isinstance(xyz, (list, tuple)) and isinstance(xyz[0], int) and len(xyz) == 3:
            x = 1
            for xs in xyz:
                x *= xs
            i = 0
            for a in ["a", "b", "c"]:
                newa = [x*xyz[i] for x in self.cell_parameter[a]]
                self.cell_parameter[a] = newa
                i += 1
                newposdir = dict()
            for atoms in self.position:
                for (i, j, k) in [(i, j, k)for i in range(0, xyz[0]) for j in range(0, xyz[1]) for k in range(0, xyz[2])]:
                    ijk = [i, j, k]
                    pos = self.position[atoms]
                    newpos = [(pos[psi]+ijk[psi])/xyz[psi]
                              for psi in range(0, 3)]
                    newposdir[(atoms[0], atoms[1]+(i*xyz[1]*xyz[2] +
                                                   j*xyz[2]+k)*self.atomtype[atoms[0]])] = newpos
            self.position = newposdir
            for atomname in self.atomtype:
                self.atomtype[atomname] *= x

        else:
            # print(type(xyz),type(xyz[0]),len(xyz))
            return None

    def findatom(self, keywords):
        # 用法举例：
        #	obj.findatom("Ba10")返回从零开始，序号为10的Ba原子
        #	obj.findatom("K")返回所有K原子
        #	obj.findatom("H8 K")返回两个列表，一个是从零开始，序号为8的H原子，一个是所有K原子
        #	obj.findatom("Ba 8")返回所有Ba原子，第二个关键字8将被忽略
        keywordslist = keywords.split()
        outlist = ["0"]
        outlistindex = 0
        outnumber = 0
        for keyword in keywordslist:
            atomname = re.search("[a-zA-Z]+", keyword)
            index = re.search("[0-9]+\-[0-9]+|[0-9]+", keyword)
            try:
                atomname = atomname.group()
            except:
                atomname = ""
            try:
                index = index.group()
            except:
                index = ""
            if atomname == "":
                continue
            elif index == "":
                start = 0
                try:
                    end = self.atomtype[atomname]-1
                except:
                    print("no such atom as %s" % atomname)
                    continue
            elif len(index.split('-')) >= 2:
                start = int(index.split('-')[0])
                end = int(index.split('-')[1])
            else:
                end = int(index)
                start = int(index)
            outnumber += 1
            if outnumber == 2:
                outlist = [outlist, ["0"]]
            elif outnumber > 2:
                outlist.append(["0"])
            outlistindex = 0
            for i in range(start, end+1):
                key = (atomname, i)
                if key in self.position:
                    if outnumber == 1:
                        try:
                            outlist[outlistindex] = key
                        except:
                            outlist.append(key)
                        outlistindex += 1
                    else:
                        try:
                            outlist[outnumber-1][outlistindex] = key
                        except:
                            outlist[outnumber-1].append(key)
                        outlistindex += 1
                else:
                    print("no such atom as %s%d" % (atomname, i))
        return outlist

    def posvector(self, atom):
        try:
            return np.array(self.position[atom])
        except:
            print("no such atom as %s!" % atom)
            return np.array([-1.0, -1.0, -1.0])

    def readqe(self, fileobject):
        strline = fileobject.readline()
        while strline[0:6] != "CELL_P":
            strline = fileobject.readline()
        for vector in ["a", "b", "c"]:
            strline = fileobject.readline()
            i = 0
            self.cell_parameter[vector] = [0.0, 0.0, 0.0]
            for num_str in strline.split():
                self.cell_parameter[vector][i] = float(num_str)
                i += 1
        # 读取晶胞参数
        strline = fileobject.readline()
        while strline[0:8] != "ATOMIC_P":
            strline = fileobject.readline()
        posmodl = strline.split('(')
        posmods = posmodl[1]
        # 提取出原子坐标类型
        self.atomtype.clear()
        self.atomtypelist.clear()
        self.position.clear()
        while True:
            # 一直读，直到没有可识别的原子位置信息行
            strline = fileobject.readline()
            postr = strline.split()
            # 把原子行分隔
            if (len(postr) != 4):
                break
            else:
                # 下面开始记录原子位置信息
                atomname = postr[0]
                # 原子的名字就是第一个
                if atomname in self.atomtype:
                    # 如果这种原子已经出现的话，只需要将它的数目加一
                    self.atomtype[atomname] += 1
                else:
                    # 如果这种原子未出现过，需要在原子类型列表中添加它，同时需要在原子类型字典中添加它，原子个数暂计为1
                    self.atomtypelist.append(atomname)
                    self.atomtype[atomname] = 1
                atomindex = self.atomtype[atomname]-1
                # 把一种原子内的编号计为atomindex，以0开始
                self.position[(atomname, atomindex)] = [0.0, 0.0, 0.0]
                # 在position字典中添加key，值待后续读取后添加
                for i in range(1, 4):
                    self.position[(atomname, atomindex)][i-1] = float(postr[i])
                    # 设置原子的位置值
        if posmods[0] == 'c':
            # c=crystal 晶体坐标
            pass
        else:
            # 如果不是晶体坐标的话，就对位置值进行转化
            array_cell_parameter = np.array(
                [self.cell_parameter["a"], self.cell_parameter["b"], self.cell_parameter["c"]])
        # 将晶胞参数转化为矩阵
            for atoms in self.position:
                # 下面开始对每一个原子的坐标进行转化
                array_xyz = np.array(self.position[atoms])
                # 将位置列表转化为矩阵
                array_xyz_out = np.dot(
                    array_xyz, np.linalg.inv(array_cell_parameter))
                # 计算
                self.position[atoms] = array_xyz_out.tolist()
                # 转回

    def readvasp(self, fileobject):
        strline = fileobject.readline()
        self.name = strline.rstrip('\n')
        # 第一行作为名称
        strline = fileobject.readline()
        # 第二行无用，不管
        self.cell_parameter.clear()
        for vector in ["a", "b", "c"]:
            strline = fileobject.readline()
            i = 0
            self.cell_parameter[vector] = [0.0, 0.0, 0.0]
            for num_str in strline.split():
                self.cell_parameter[vector][i] = float(num_str)
                i += 1
        # 读取晶胞参数
        strline = fileobject.readline()
        atomnames = strline
        strline = fileobject.readline()
        atomnumbers = strline
        numberlist = atomnumbers.split()
        self.atomtype.clear()
        self.position.clear()
        i = 0
        self.atomtypelist = atomnames.split()
        types = len(self.atomtypelist)
        # print(self.atomtypelist)
        # print(numberlist)
        for atomname in self.atomtypelist:
            # print("%d-%d"%(i,int(numberlist[i])))
            self.atomtype[atomname] = int(numberlist[i])
            # self.atomtypelist.append(atomname)
            i += 1
            if i >= types:
                break
        # 读取原子个数
        strline = fileobject.readline()
        # 读取模式，D-笛卡尔坐标or晶体坐标，C-直接坐标
        self.position.clear()
        if strline[0] == "D":
            for atomname in atomnames.split():
                i = 0
                number = self.atomtype[atomname]
                while i < number:
                    strline = fileobject.readline()
                    j = 0
                    self.position[(atomname, i)] = [0.0, 0.0, 0.0]
                    for num_str in strline.split():
                        self.position[(atomname, i)][j] = float(num_str)
                        j += 1
                    i += 1
        elif strline[0] == "C":
            array_cell_parameter = np.array(
                [self.cell_parameter["a"], self.cell_parameter["b"], self.cell_parameter["c"]])
            for atomname in atomnames.splict():
                i = 0
                number = self.atomtype[atomname]
                while i < number:
                    strline = fileobject.readline()
                    j = 0
                    xyz = [0.0, 0.0, 0.0]
                    for num_str in strline.split():
                        xyz[j] = float(num_str)
                        j += 1
                    array_xyz = np.array(xyz)
                    array_xyz_out = dot(array_xyz, inv(array_cell_parameter))
                    self.position[(atomname, i)] = array_xyz_out.tolist()
                    i += 1

    def replace(self, atomSrc, strDis):
        # 在位置字典中中搜索atomSrc，将改原子改名为strdis
        strDis = strDis.strip()
        # if isinstance(atomSrc,list):
        #	out=self.copy()
        #	for atoms in atomSrc:
        #		out=out.replace(atoms,strDis)
        #	return out
        if atomSrc in self.position:
            out = self.copy()
            # print("_______________________")
            # print(out)
            if strDis in out.atomtypelist:
                # 目标原子如果已存在的话，其总数加一
                index = out.atomtype[strDis]
                out.atomtype[strDis] += 1
            else:
                # 目标原子不存在的话，向列表中添加它
                out.atomtype[strDis] = 1
                out.atomtypelist.append(strDis)
                index = 0
            out.position[(strDis, index)] = out.position[atomSrc]
            delname = atomSrc[0]
            del out.position[atomSrc]
            j = 0
            for i in range(0, out.atomtype[delname]+1):
                try:
                    poslist = out.position[(delname, i)]
                    del out.position[(delname, i)]
                    out.position[(delname, j)] = poslist
                    j += 1
                except:
                    pass
            if out.atomtype[delname] == 1:
                out.atomtypelist.remove(delname)
                del out.atomtype[delname]
            else:
                out.atomtype[delname] -= 1
            return out
        else:
            return None

    def writevasp(self, fileobject):
        fileobject.write("%s\n" % self.name)
        fileobject.write("1.0\n")
        for a in ["a", "b", "c"]:
            fileobject.write("\t%s\t%s\t%s\n" % tuple(self.cell_parameter[a]))
        for key in self.atomtypelist:
            fileobject.write("\t%s" % key)
        fileobject.write("\n")
        for key in self.atomtypelist:
            fileobject.write("\t%d" % self.atomtype[key])
        fileobject.write("\nDirect\n")
        for key in self.atomtypelist:
            for index in range(0, self.atomtype[key]):
                fileobject.write("\t%s\t%s\t%s\n" %
                                 tuple(self.position[(key, index)]))

    def writeqe(self, fileobject):
        fileobject.write("CELL_PARAMETERS (angstrom)\n")
        for v in ["a", "b", "c"]:
            outstr = ""
            paralist = self.cell_parameter[v]
            for i in [0, 1, 2]:
                outstr = outstr+"\t%s" % paralist[i]
            fileobject.write(outstr+"\n")
        fileobject.write("ATOMIC_POSITIONS (crystal)\n")
        for atomtype in self.atomtypelist:
            for i in range(0, self.atomtype[atomtype]):
                atoms = (atomtype, i)
                outstr = atomtype
                for j in [0, 1, 2]:
                    outstr = outstr+"\t%s" % (self.position[atoms][j])
                fileobject.write(outstr+"\n")


def testit():
    testcrystal = crystal()
    print(testcrystal)


def help():
    helppage = '''
	crystal:
		晶胞类,包括一下几个基本属性
	name：字符串，定义晶胞名称
	cell_parameter: 字典，定义晶胞参数，key为"a"、"b"、"c"，value为列表
	atomtype: 字典，定义晶胞中原子总类和每种的数目，key为原子名，value为原子个数
	atomtypelist: 列表，定义晶胞中原子的种类列表
	position: 字典，定义晶胞中各原子位置，key为元组，由原子名+序号构成，下称“晶胞原子组”，序号以0开始
	基本方法：
	atomdistance(self,left,right):
		left，right是元组，如果两者任何一个不在原子位置字典position中，返回-1，否则返回两个原子间的距离
	CellParameter(self):
		返回以numpy中矩阵形式存储的晶胞参数
	check(self):
		检查晶胞是否“自洽”，正常返回1,不正常返回一个小于等于零的数：
		0：原子总数与位置字典中key数不匹配
		……
	clear(self):
		将所有的属性置空
	findatom(self,keywords):
		根据keywords在晶胞中查找原子，keywords是字符串，返回一个列表。
		如果只有一个“有效”的keyword，则列表的元素是符合条件的“晶胞原子组”。
		如果有多个“有效”的keywords，则列表的元素是每个“有效”keywords匹配得到的列表。
		合法的keyword形式有三种： K Ba5 O23-40
	posvector(self,atom):
		atom是“晶胞原子组”，返回一个以numpy中矩阵形式存储的输入atom所指的位置向量
	readqe(self,fileobject):
		从qe可读格式的文件中读取晶胞
	readvasp(self,fileobject):
		从vasp格式的文件中读取晶胞
	writevasp(self,fileobject):
		将晶胞输出为vasp格式的文件
	writeqe(self,fileobject):
		将晶胞输出为qe可读格式的文件
	'''
    print(helppage)
