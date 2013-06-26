#include "pw_flow.h"

namespace PowerSystem
{
	BusAndBranch::BusAndBranch()
	{	
	}

	void BusAndBranch::FromConfigFile(string filename)
	{
		fstream fp;
		fp.open(filename.c_str(), ios::in);

		if (!fp.is_open())
		{
			throw _Exception("file open error(In FromConfigFile)");
		}

		char line[80];
		fp.getline(line, 78);
		if (string(line) != "LOVE_YP_POWERFLOW_DATA")
		{
			throw _Exception("Miss string 'LOVE_YP_POWERFLOW_DATA'as the first one in file(In FromConfigFile)");
		}

		while (!fp.eof())
		{
			fp.getline(line, 78);
			string text = ExString::Trim(string(line));
			if (text.find("BUS_INFO_BEGIN") != -1)
			{
				string del[] = {string(":")};
				bus.reserve(atoi(ExString::Splite(text, StrArray(del, del + 1))[1].c_str()));
				while (!fp.eof())
				{
					fp.getline(line, 78);
					text = ExString::Trim(string(line));
					if (text.find("BUS_INFO_END") != -1)
					{
						break;
					}
					if (text.length() == 0 || text[0] == '#')
					{
						continue;
					}else{
						string sp[] = {"\t"};
						StrArray split = ExString::Splite(text, StrArray(sp, sp + 1));
						if (split.size() != 13)
						{
							throw _Exception("Bus data in the file is not correct(In FromConfigFile)");
						}
						bus.push_back(Bus(atoi(split[0].c_str()), split[1], Bus::NodeType(atoi(split[2].c_str())), atof(split[3].c_str()), atof(split[4].c_str()), atof(split[5].c_str()), atof(split[6].c_str()), atof(split[7].c_str()), atof(split[8].c_str()), atof(split[9].c_str()), atof(split[10].c_str()), atof(split[11].c_str()), atof(split[12].c_str())));
					}
				}
			}

			if (text.find("BRANCH_INFO_BEGIN") != -1)
			{
				string del[] = {string(":")};
				branch.reserve(atoi(ExString::Splite(text, StrArray(del, del + 1))[1].c_str()));
				while (!fp.eof())
				{
					fp.getline(line, 78);
					text = ExString::Trim(string(line));
					if (text.find("BRANCH_INFO_END") != -1)
					{
						break;
					}
					if (text.length() == 0 || text[0] == '#')
					{
						continue;
					}else{
						string sp[] = {"\t"};
						StrArray split = ExString::Splite(text, StrArray(sp, sp + 1));
						if (split.size() != 8)
						{
							throw _Exception("Branch data in the file is not correct(In FromConfigFile)");
						}
						branch.push_back(Branch(atoi(split[0].c_str()), atoi(split[1].c_str()), split[2], (Branch::BranchType)atoi(split[3].c_str()), atof(split[4].c_str()), atof(split[5].c_str()), atof(split[6].c_str()), atof(split[7].c_str())));
					}
				}
			}

			if (ExString::Trim(text).length() == 0 || ExString::Trim(text)[0] == '#')
			{
				continue;
			}
		}
		//PrintBuses();
		//PrintBranches();
		//system("pause");
		fp.clear();
		fp.close();		
	}

	void BusAndBranch::PrintBranch(const Branch& one)
	{
		cout << one.i << "\t" << one.j << '\t' << one.Name << '\t' << one.Type << '\t' << one.R << '\t' << one.X << '\t' << one.B << '\t' << one.K << endl;
	}

	void BusAndBranch::PrintBus(const Bus& one)
	{
		cout << one.i << "\t" << one.Name << "\t" << (int)one.Type << "\t" << one.V << "\t" << one.cita << "\t" << one.BS << "\t"  << one.GP  << "\t" << one.GQ << "\t" << one.LP << "\t" << one.LQ  << "\t" << one.Qmax << "\t" << one.Qmin << "\t" << one.b << '\t' << one.P << "\t" << one.Q << endl;
	}

	void BusAndBranch::PrintBranches() const
	{
		for (vector<Branch>::const_iterator it = branch.begin(); it != branch.end(); it++)
		{
			PrintBranch(*it);
		}
	}

	void BusAndBranch::PrintBuses() const
	{
		for (vector<Bus>::const_iterator it = bus.begin(); it != bus.end(); it++)
		{
			PrintBus(*it);
		}
	}

	int BusAndBranch::GetBusSize() const
	{
		return bus.size();
	}

	int BusAndBranch::GetBranchSize() const
	{
		return branch.size();
	}

	BusAndBranch::Bus* BusAndBranch::GetBus(int id)
	{
		if (id >= (int)bus.size())
			throw _Exception("Index out of the range(In GetBus)");
		return &bus[id];
	}

	BusAndBranch::Branch* BusAndBranch::GetBranch(int id)
	{
		if (id >= (int)branch.size())
			throw _Exception("Index out of the range(In GetBranch)");
		return &branch[id];
	}

	_PowerFlow::_PowerFlow()
	{
		nPQ = 0;
		nPV = 0;
		nSlack = 0;
	}

	_PowerFlow::~_PowerFlow()
	{
		Y.clear();
		NodeIdMap.clear();
	}

	void _PowerFlow::InitBABFromFile(string filename)
	{
		BAB.FromConfigFile(filename);
		Y.reserve(BAB.GetBusSize() + 2 * BAB.GetBranchSize());
		NodeIdMap.reserve(BAB.GetBusSize());
	}

	int _PowerFlow::GetNodeNum() const
	{
		return NodeIdMap.size();
	}

	void _PowerFlow::SortBABByType()
	{
		//No efficiency, move too much memory
		/*for (int i = 0; i<BAB.GetBusSize() - 1; i++)
		{
			for (int j = i + 1; j<BAB.GetBusSize(); j++)
			{
				if (int(BAB.GetBus(i)->Type) > int(BAB.GetBus(j)->Type))
				{
					BusAndBranch::Bus temp;
					temp = *BAB.GetBus(i);
					*BAB.GetBus(i) = *BAB.GetBus(j);
					*BAB.GetBus(j) = temp;
				}
			}
		}*/

		for (int i = 0; i<BAB.GetBusSize(); i++)
		{
			NodeIdMap.push_back(Pair(int(BAB.GetBus(i)->Type), i));

			switch(BAB.GetBus(i)->Type)
			{
			case BusAndBranch::Bus::SLACK:
				nSlack++;
				break;
			case BusAndBranch::Bus::PQ:
				nPQ++;
				break;
			case BusAndBranch::Bus::PV:
				nPV++;
				break;
			default:
				throw _Exception("Unknown bus type(In SortBABByType)");
				break;
			}
		}

		if (nSlack != 1)
			throw _Exception("The number of slack buses is not '1'(In SortBABByType)");

		for (vector<Pair>::iterator it = NodeIdMap.begin(); it < NodeIdMap.end() - 1; it++)
		{
			for (vector<Pair>::iterator itj = it + 1; itj != NodeIdMap.end(); itj++)
			{
				if (it->i * NodeIdMap.size() + it->j > itj->i * NodeIdMap.size() + itj->j )
				{
					Pair temp = *it;
					*it = *itj;
					*itj = temp;
				}
			}
		}

		for (vector<Pair>::iterator it = NodeIdMap.begin(); it != NodeIdMap.end(); it++)
		{
			it->i = distance(NodeIdMap.begin(), it);
		}
	}

	int _PowerFlow::GetSystemMap(int i) const
	{
		if (i >= (int)NodeIdMap.size())
		{
			throw _Exception("Index out of the NodeIdMap's range(In GetSystemMap)");
		}

		return NodeIdMap[i].j;
	}

	int _PowerFlow::GetVirtualId(int i)
	{
		for (vector<Pair>::const_iterator it = NodeIdMap.begin(); it != NodeIdMap.end(); it++)
		{
			if (BAB.GetBus(it->j)->i == i)
			{
				return it->i;
			}
		}
		return -1;
	}

	double _PowerFlow::GetPi(int i, const vector<double>& V, const vector<double>& cita) const
	{
		double ret = 0.0;

		for (vector<ThreeElements<complex<double>>>::const_iterator it = Y.begin(); it != Y.end(); it++)
		{
			if (it->i == i)
			{
				int j = it->j;				
				ret += V[j] * (it->data.real() * cos (cita[i] - cita[j]) + it->data.imag() * sin(cita[i] - cita[j]));				
			}
		}
		return ret * V[i];
	}

	double _PowerFlow::GetQi(int i, const vector<double>& V, const vector<double>& cita) const
	{
		double ret = 0.0;

		for (vector<ThreeElements<complex<double>>>::const_iterator it = Y.begin(); it != Y.end(); it++)
		{
			if (it->i == i)
			{
				int j = it->j;	
				ret += V[j] * (it->data.real() * sin(cita[i] - cita[j]) - it->data.imag() * cos (cita[i] - cita[j]));
			}
		}
		return ret * V[i];
	}

	void _PowerFlow::BranchFlowCal()
	{
		Pij.reserve(BAB.GetBranchSize() * 2);
		Qij.reserve(BAB.GetBranchSize() * 2);
		for (int i = 0; i<BAB.GetBranchSize(); i++)
		{
			int _i = BAB.GetBranch(i)->i, _j = BAB.GetBranch(i)->j;
			int __i = GetVirtualId(_i), __j = GetVirtualId(_j);
		
			BusAndBranch::Bus *pBuI = BAB.GetBus(GetSystemMap(__i)), *pBuJ = BAB.GetBus(GetSystemMap(__j));
			double gij = 0.0, bij = 0.0;
			for (vector<ThreeElements<complex<double>>>::const_iterator it = Y.begin(); it != Y.end(); it++)
			{
				if (it->i == __i && it->j == __j)
				{
					gij = -it->data.real();
					bij = -it->data.imag();
					break;
				}

			}

			double g0i = 0, g0j = 0;
			for (int k = 0; k < BAB.GetBranchSize(); k++)
			{
				BusAndBranch::Branch* pBr = BAB.GetBranch(k);
				if (pBr->j == _j  && pBr->i == _i)
				{
					if (pBr->Type == BusAndBranch::Branch::TRANSFORMER)
					{
						g0i = gij * (pBr->K - 1);
						g0j = gij * (1 - pBr->K) / pBr->K;
						break;
					}		
				}
			}

			Pij.push_back(ThreeElements<double>(_i, _j, pow(pBuI->V, 2.0) * (gij + g0i)- pBuI->V * pBuJ->V * (gij * cos(pBuI->cita - pBuJ->cita) + bij * sin(pBuI->cita - pBuJ->cita))));
			Pij.push_back(ThreeElements<double>(_j, _i, pow(pBuJ->V, 2.0) * (gij + g0j) - pBuI->V * pBuJ->V * (gij * cos(pBuI->cita - pBuJ->cita) - bij * sin(pBuI->cita - pBuJ->cita))));
		
			double b0i = 0.0, b0j = 0;;
			for (int k = 0; k < BAB.GetBranchSize(); k++)
			{
				BusAndBranch::Branch* pBr = BAB.GetBranch(k);
				if (pBr->i == _i  && pBr->j == _j)
				{
					if (pBr->Type == BusAndBranch::Branch::TRANSFORMER)
					{
						b0i = bij * (pBr->K - 1);
						b0j = bij * (1 - pBr->K) / pBr->K;
					}else{
						b0i = b0j = pBr->B;
					}
					break;
				}
			}

			Qij.push_back(ThreeElements<double>(_i, _j, -1 * pow(pBuI->V, 2.0) * (bij + b0i) + pBuI->V * pBuJ->V * (bij * cos(pBuI->cita - pBuJ->cita) - gij * sin(pBuI->cita - pBuJ->cita))));	
			Qij.push_back(ThreeElements<double>(_j, _i, -1 * pow(pBuJ->V, 2.0) * (bij + b0j) + pBuI->V * pBuJ->V * (bij * cos(pBuI->cita - pBuJ->cita) + gij * sin(pBuI->cita - pBuJ->cita))));
		}
	}

	FastPowerFlow::FastPowerFlow()
	{
	}

	FastPowerFlow::~FastPowerFlow()
	{
	}

	void FastPowerFlow::CreateYMatrix()
	{
		SortBABByType();
		for (vector<Pair>::const_iterator it = NodeIdMap.begin(); it != NodeIdMap.end(); it++)
		{
			int realId = BAB.GetBus(it->j)->i;
			Y.push_back(ThreeElements<complex<double>>(it->i, it->i, complex<double>(0, BAB.GetBus(it->j)->b)));
			int major = Y.size() - 1;

			for (int i = 0; i < BAB.GetBranchSize(); i++)
			{
				BusAndBranch::Branch* pBr = BAB.GetBranch(i);

				if (pBr->i == realId || pBr->j == realId)
				{
					bool left = pBr->i == realId;
					complex<double> y = complex<double>(-1.0) / complex<double>(pBr->R , pBr->X), ym = complex<double>(-1.0) * y;
					switch (pBr->Type)
					{
					case BusAndBranch::Branch::LINE:
						break;
					case BusAndBranch::Branch::TRANSFORMER:
						y /= complex<double>(pBr->K);
						if (!left)
						{
							ym /= complex<double>(pow(pBr->K, 2.0));
						}
						break;
					default:
						throw _Exception("Unknown branch type(In CreateYMatrix)");
						break;
					}
					if (GetVirtualId(pBr->i) >= it->i && GetVirtualId(pBr->j) >= it->i)
					{
						Y.push_back(ThreeElements<complex<double>>(GetVirtualId(pBr->i), GetVirtualId(pBr->j), y));
						Y.push_back(ThreeElements<complex<double>>(GetVirtualId(pBr->j), GetVirtualId(pBr->i), y));
					}
					
					Y[major].data += ym + complex<double>(0.0, pBr->B);

				}
			}
		}
		SparseMatrix<complex<double>>::SortTEs(Y, GetNodeNum());
		SparseMatrix<complex<double>> sm(Y, Y.size(), GetNodeNum());
		Matrix<complex<double>> m;
		sm.ToMatrix(m);
		m.ToFile("YY.txt");		
	}

	void FastPowerFlow::CreateB1_XB()
	{
		vector<ThreeElements<double>> _B1;
		int offset = nSlack;
		_B1.reserve(BAB.GetBusSize() + BAB.GetBranchSize() * 2);
		B1Map.reserve(nPQ + nPV);

		vector<int> index;
		index.reserve(nPQ + nPV);
		for (vector<ThreeElements<complex<double>>>::const_iterator it = Y.begin(); it != Y.end(); it++)
		{
			int i = it->i - offset;
			int j = it->j - offset;
			if (i >= 0 && j >= 0)
			{
				_B1.push_back(ThreeElements<double>(i, j, 1.0 / (complex<double>(1) / it->data).imag()));
				if (i == j)
				{
					B1Map.push_back(it->i);
					index.push_back(_B1.size() - 1);
				}
			}
		}

		for (vector<int>::const_iterator it = index.begin(); it != index.end(); it++)
		{
			int realId = BAB.GetBus(GetSystemMap(_B1[*it].i + offset))->i;
			_B1[*it].data = 0;
			for (int i = 0; i < BAB.GetBranchSize(); i++)
			{
				BusAndBranch::Branch* pBr = BAB.GetBranch(i);
				
				if (realId == pBr->i || realId == pBr->j)
				{
					_B1[*it].data += 1.0 / pBr->X;
				}
			}
		}

		int add = 0;
		if (nPQ + nPV > 2)
		{
			Graph<double> gh;
			gh.FromTEs(_B1);
			vector<int> ps;		
			vector<PowerSystem::Pair> psAdd;	
			gh.Optimization(ps, add, &psAdd);
			add *= 2;
			gh.Print();
			PowerSystem::test_Print(ps);
			cout << "\n注入元个数" << add<< endl;
			test_Print(psAdd);
			vector<int> temp;
			for (vector<int>::iterator it = ps.begin(); it != ps.end(); it++)
			{
				temp.push_back(B1Map[*it]);
			}
			B1Map.assign(temp.begin(), temp.end());

			for (vector<ThreeElements<double>>::iterator it = _B1.begin(); it !=_B1.end(); it++)
			{
				bool s1 = false, s2 = false;

				for (int i = 0; i<(int)ps.size(); i++)
				{
					if (!s1 && it->i == ps[i])
					{
						it->i = i;
						s1 = true;
					}

					if (!s2 && it->j == ps[i])
					{
						it->j = i;
						s2 = true;
					}

					if (s1 && s2)
					{
						break;
					}
				}
			}
			SparseMatrix<double>::SortTEs(_B1, nPQ + nPV);
		}

		B1.FromTEs(_B1, _B1.size() + add, nPQ + nPV);
		//Matrix<double> m;
		//B1.ToMatrix(m);
		//m.Print();

	}

	void FastPowerFlow::CreateB2_XB()
	{
		vector<ThreeElements<double>> _B2;
		int offset = nSlack + nPV;
		_B2.reserve(BAB.GetBusSize() + BAB.GetBranchSize() * 2);
		B2Map.reserve(nPQ);

		for (vector<ThreeElements<complex<double>>>::const_iterator it = Y.begin(); it != Y.end(); it++)
		{
			int i = it->i - offset;
			int j = it->j - offset;
			if (i >= 0 && j >= 0)
			{
				_B2.push_back(ThreeElements<double>(i, j, -1 * it->data.imag()));
				if (i == j)
				{
					B2Map.push_back(it->i);
				}
			}
		}

		//test_Print(_B2);
		int add = 0;
		if (nPQ > 2)
		{
			Graph<double> gh;
			gh.FromTEs(_B2);
			vector<int> ps;		
			vector<PowerSystem::Pair> psAdd;	
			gh.Optimization(ps, add, &psAdd);
			add *= 2;
			//gh.Print();
			//PowerSystem::test_Print(ps);
			//cout << "\n注入元个数" << add * 2<< endl;
			//test_Print(psAdd);
			vector<int> temp;
			for (vector<int>::iterator it = ps.begin(); it != ps.end(); it++)
			{
				temp.push_back(B2Map[*it]);
			}
			B2Map.assign(temp.begin(), temp.end());

			for (vector<ThreeElements<double>>::iterator it = _B2.begin(); it !=_B2.end(); it++)
			{
				bool s1 = false, s2 = false;

				for (int i = 0; i<(int)ps.size(); i++)
				{
					if (!s1 && it->i == ps[i])
					{
						it->i = i;
						s1 = true;
					}

					if (!s2 && it->j == ps[i])
					{
						it->j = i;
						s2 = true;
					}

					if (s1 && s2)
					{
						break;
					}
				}
			}
			//test_Print(_B2);
			SparseMatrix<double>::SortTEs(_B2, nPQ);
		}

		//test_Print(_B2);
		B2.FromTEs(_B2, _B2.size() + add, nPQ);
		//Matrix<double> m;
		//B2.ToMatrix(m);
		//m.Print();
	}

	void FastPowerFlow::CreateB1_BX()
	{
		vector<ThreeElements<double>> _B1;
		int offset = nSlack;
		_B1.reserve(BAB.GetBusSize() + BAB.GetBranchSize() * 2);
		B1Map.reserve(nPQ + nPV);

		vector<int> index;
		index.reserve(nPQ + nPV);
		for (vector<ThreeElements<complex<double>>>::const_iterator it = Y.begin(); it != Y.end(); it++)
		{
			int i = it->i - offset;
			int j = it->j - offset;
			if (i >= 0 && j >= 0)
			{
				_B1.push_back(ThreeElements<double>(i, j, -1 * it->data.imag()));
				for (int b = 0; b < BAB.GetBranchSize(); b++)
				{
					BusAndBranch::Branch *pBr = BAB.GetBranch(b);
					int virtualI = GetVirtualId(pBr->i), virtualJ = GetVirtualId(pBr->j);
					if ((virtualI == it->i && virtualJ == it->j || virtualI == it->j && virtualJ == it->i) && pBr->Type == BusAndBranch::Branch::TRANSFORMER)
					{
						_B1[_B1.size() - 1].data *= pBr->K;
						break;
					}
				}
				if (i == j)
				{
					B1Map.push_back(it->i);
					index.push_back(_B1.size() - 1);
				}
			}
		}

		for (vector<int>::const_iterator it = index.begin(); it != index.end(); it++)
		{
			_B1[*it].data = 0;
			for (vector<ThreeElements<complex<double>>>::const_iterator pY = Y.begin(); pY != Y.end(); pY++)
			{
				if (pY->i == _B1[*it].i + offset && pY->j != _B1[*it].j + offset)
				{		
					double data = pY->data.imag();
					for (int b = 0; b < BAB.GetBranchSize(); b++)
					{
						BusAndBranch::Branch *pBr = BAB.GetBranch(b);
						int virtualI = GetVirtualId(pBr->i), virtualJ = GetVirtualId(pBr->j);
						if ((virtualI == pY->i && virtualJ == pY->j || virtualI == pY->j && virtualJ == pY->i) && pBr->Type == BusAndBranch::Branch::TRANSFORMER)
						{
							data *= pBr->K;
							break;
						}
					}
					_B1[*it].data += data;
				}
			}
		}

		int add = 0;
		if (nPQ + nPV > 2)
		{
			Graph<double> gh;
			gh.FromTEs(_B1);
			vector<int> ps;		
			vector<PowerSystem::Pair> psAdd;	
			gh.Optimization(ps, add, &psAdd);
			add *= 2;
			vector<int> temp;
			for (vector<int>::iterator it = ps.begin(); it != ps.end(); it++)
			{
				temp.push_back(B1Map[*it]);
			}
			B1Map.assign(temp.begin(), temp.end());

			for (vector<ThreeElements<double>>::iterator it = _B1.begin(); it !=_B1.end(); it++)
			{
				bool s1 = false, s2 = false;

				for (int i = 0; i<(int)ps.size(); i++)
				{
					if (!s1 && it->i == ps[i])
					{
						it->i = i;
						s1 = true;
					}

					if (!s2 && it->j == ps[i])
					{
						it->j = i;
						s2 = true;
					}

					if (s1 && s2)
					{
						break;
					}
				}
			}

			SparseMatrix<double>::SortTEs(_B1, nPQ + nPV);
		}

		B1.FromTEs(_B1, _B1.size() + add, nPQ + nPV);
		Matrix<double> m;
		B1.ToMatrix(m);
		m.ToFile();
	}

	void FastPowerFlow::CreateB2_BX()
	{
		vector<ThreeElements<double>> _B2;
		int offset = nSlack + nPV;
		_B2.reserve(BAB.GetBusSize() + BAB.GetBranchSize() * 2);
		B2Map.reserve(nPQ);

		vector<int> index;
		index.reserve(nPQ);
		for (vector<ThreeElements<complex<double>>>::const_iterator it = Y.begin(); it != Y.end(); it++)
		{
			int i = it->i - offset;
			int j = it->j - offset;
			if (i >= 0 && j >= 0)
			{
				//_B2.push_back(ThreeElements<double>(i, j, 1.0 / (complex<double>(1) / it->data).imag()));
				for (int b = 0; b < BAB.GetBranchSize(); b++)
				{
					BusAndBranch::Branch *pBr = BAB.GetBranch(b);
					int virtualI = GetVirtualId(pBr->i), virtualJ = GetVirtualId(pBr->j);
					if (virtualI == it->i && virtualJ == it->j || virtualI == it->j && virtualJ == it->i)
					{
						_B2.push_back(ThreeElements<double>(i, j, - 1.0 / pBr->X));
						break;
					}
				}
				if (i == j)
				{
					B2Map.push_back(it->i);
					_B2.push_back(ThreeElements<double>(i, j, 0));
					index.push_back(_B2.size() - 1);
				}
			}
		}

		for (vector<int>::const_iterator it = index.begin(); it != index.end(); it++)
		{ 
			BusAndBranch::Bus *pBu = BAB.GetBus(GetSystemMap(_B2[*it].i + offset));
			int realId = pBu->i;
			_B2[*it].data = -pBu->b;
			for (int i = 0; i < BAB.GetBranchSize(); i++)
			{
				BusAndBranch::Branch* pBr = BAB.GetBranch(i);
				
				if (realId == pBr->i || realId == pBr->j)
				{
					if (pBr->Type == BusAndBranch::Branch::TRANSFORMER && realId == pBr->j)
					{
						_B2[*it].data += 1.0 / pBr->X / pow(pBr->K, 2.0) - pBr->B;
					}else{
						_B2[*it].data += 1.0 / pBr->X - pBr->B;
					}
				}
			}
		}


		int add = 0;
		if (nPQ > 2)
		{
			Graph<double> gh;
			gh.FromTEs(_B2);
			vector<int> ps;		
			vector<PowerSystem::Pair> psAdd;	
			gh.Optimization(ps, add, &psAdd);
			add *= 2;
			vector<int> temp;
			for (vector<int>::iterator it = ps.begin(); it != ps.end(); it++)
			{
				temp.push_back(B2Map[*it]);
			}
			B2Map.assign(temp.begin(), temp.end());

			for (vector<ThreeElements<double>>::iterator it = _B2.begin(); it !=_B2.end(); it++)
			{
				bool s1 = false, s2 = false;

				for (int i = 0; i<(int)ps.size(); i++)
				{
					if (!s1 && it->i == ps[i])
					{
						it->i = i;
						s1 = true;
					}

					if (!s2 && it->j == ps[i])
					{
						it->j = i;
						s2 = true;
					}

					if (s1 && s2)
					{
						break;
					}
				}
			}
			SparseMatrix<double>::SortTEs(_B2, nPQ);
		}

		B2.FromTEs(_B2, _B2.size() + add, nPQ);
		Matrix<double> m;
		B2.ToMatrix(m);
		m.ToFile();
	}

	void FastPowerFlow::CreateB1_PQ()
	{
		vector<ThreeElements<double>> _B1;
		int offset = nSlack;
		_B1.reserve(BAB.GetBusSize() + BAB.GetBranchSize() * 2);
		B1Map.reserve(nPQ + nPV);

		for (vector<ThreeElements<complex<double>>>::const_iterator it = Y.begin(); it != Y.end(); it++)
		{
			int i = it->i - offset;
			int j = it->j - offset;
			if (i >= 0 && j >= 0)
			{
				_B1.push_back(ThreeElements<double>(i, j, -1 * it->data.imag()));
				if (i == j)
				{
					B1Map.push_back(it->i);
				}
			}
		}

		int add = 0;
		if (nPQ + nPV > 2)
		{
			Graph<double> gh;
			gh.FromTEs(_B1);
			vector<int> ps;		
			vector<PowerSystem::Pair> psAdd;	
			gh.Optimization(ps, add, &psAdd);
			add *= 2;
			vector<int> temp;
			for (vector<int>::iterator it = ps.begin(); it != ps.end(); it++)
			{
				temp.push_back(B1Map[*it]);
			}
			B1Map.assign(temp.begin(), temp.end());

			for (vector<ThreeElements<double>>::iterator it = _B1.begin(); it !=_B1.end(); it++)
			{
				bool s1 = false, s2 = false;

				for (int i = 0; i<(int)ps.size(); i++)
				{
					if (!s1 && it->i == ps[i])
					{
						it->i = i;
						s1 = true;
					}

					if (!s2 && it->j == ps[i])
					{
						it->j = i;
						s2 = true;
					}

					if (s1 && s2)
					{
						break;
					}
				}
			}

			SparseMatrix<double>::SortTEs(_B1, nPQ + nPV);
		}

		B1.FromTEs(_B1, _B1.size() + add, nPQ + nPV);
	}

	void FastPowerFlow::CreateB2_PQ()
	{
		vector<ThreeElements<double>> _B2;
		int offset = nSlack + nPV;
		_B2.reserve(BAB.GetBusSize() + BAB.GetBranchSize() * 2);
		B2Map.reserve(nPQ);

		for (vector<ThreeElements<complex<double>>>::const_iterator it = Y.begin(); it != Y.end(); it++)
		{
			int i = it->i - offset;
			int j = it->j - offset;
			if (i >= 0 && j >= 0)
			{
				_B2.push_back(ThreeElements<double>(i, j, -1 * it->data.imag()));
				if (i == j)
				{
					B2Map.push_back(it->i);
				}
			}
		}

		int add = 0;
		if (nPQ > 2)
		{
			Graph<double> gh;
			gh.FromTEs(_B2);
			vector<int> ps;		
			vector<PowerSystem::Pair> psAdd;	
			gh.Optimization(ps, add, &psAdd);
			add *= 2;
			vector<int> temp;
			for (vector<int>::iterator it = ps.begin(); it != ps.end(); it++)
			{
				temp.push_back(B2Map[*it]);
			}
			B2Map.assign(temp.begin(), temp.end());

			for (vector<ThreeElements<double>>::iterator it = _B2.begin(); it !=_B2.end(); it++)
			{
				bool s1 = false, s2 = false;

				for (int i = 0; i<(int)ps.size(); i++)
				{
					if (!s1 && it->i == ps[i])
					{
						it->i = i;
						s1 = true;
					}

					if (!s2 && it->j == ps[i])
					{
						it->j = i;
						s2 = true;
					}

					if (s1 && s2)
					{
						break;
					}
				}
			}
			SparseMatrix<double>::SortTEs(_B2, nPQ);
		}

		B2.FromTEs(_B2, _B2.size() + add, nPQ);
		B2.FromTEs(_B2, _B2.size() + add, nPQ);
	}

	void FastPowerFlow::GetB1FactorTable()
	{
		B1.LDUDecomposition(false);
	}

	void FastPowerFlow::GetB2FactorTable()
	{
		B2.LDUDecomposition(false);
	}

	void FastPowerFlow::CreateB1()
	{
		switch (type)
		{
		case XB:
			CreateB1_XB();
			break;
		case BX:
			CreateB1_BX();
			break;
		case PQ:
			CreateB1_PQ();
			break;
		default:
			throw _Exception("No this kind of algorithm(In CreateB1)");
			break;
		}
	}

	void FastPowerFlow::CreateB2()
	{
		switch (type)
		{
		case XB:
			CreateB2_XB();
			break;
		case BX:
			CreateB2_BX();
			break;
		case PQ:
			CreateB2_PQ();
			break;
		default:
			throw _Exception("No this kind of algorithm(In CreateB2)");
			break;
		}
	}


	void FastPowerFlow::FlowCal(int type, int times, double e)
	{	
		B1.Clear();
		B2.Clear();
		B1Map.clear();
		B2Map.clear();
		Pij.clear();
		Qij.clear();
		this->type = (Type)type;	
		CreateB1();
		CreateB2();
		GetB1FactorTable();
		GetB2FactorTable();
		vector<double> dP(nPQ + nPV), P(nPQ + nPV), cita(nPQ + nPV, 0.0), cV(nPV, 0), V(nPQ, 1.0), dQ(nPQ), Q(nPQ), PV(nPQ + nPV);
		vector<double> sV(GetNodeNum()), sCita(GetNodeNum(), 0.0);
		sV[0] = BAB.GetBus(GetSystemMap(0))->V;
		sCita[0] = BAB.GetBus(GetSystemMap(0))->cita;

		int i = 0;
		for (vector<int>::const_iterator it = B1Map.begin(); it != B1Map.end(); it++)
		{
			BusAndBranch::Bus *pBu = BAB.GetBus(GetSystemMap(*it));
			P[distance<vector<int>::const_iterator>(B1Map.begin(), it)] = pBu->GP - pBu->LP;
		}

		for (vector<int>::const_iterator it = B2Map.begin(); it != B2Map.end(); it++)
		{
			BusAndBranch::Bus *pBu = BAB.GetBus(GetSystemMap(*it));
			Q[distance<vector<int>::const_iterator>(B2Map.begin(), it)] = pBu->GQ - pBu->LQ;
		}
	
		int k = 0;
		int KP = 1, KQ = 1;
		cout << k << endl;
		while(k <= times)
		{
			for (vector<int>::const_iterator it = B1Map.begin(); it != B1Map.end(); it++)
			{	
				int _i = distance<vector<int>::const_iterator>(B1Map.begin(), it);
				BusAndBranch::Bus *pBu = BAB.GetBus(GetSystemMap(*it));
				if (pBu->Type == BusAndBranch::Bus::PV)
				{
					PV[_i] = pBu->V;
				}else{
					PV[_i] = V[GetB2Virtual(*it)];	
				}
				sV[*it] = PV[_i];
				//cout << "V" << *it << ':' << sV[*it] << endl ;
			}
			
			double maxPe = -1;
			for (vector<int>::const_iterator it = B1Map.begin(); it != B1Map.end(); it++)
			{	
				int _i = distance<vector<int>::const_iterator>(B1Map.begin(), it);
				dP[_i] = (P[_i] - GetPi(*it, sV, sCita)) / PV[_i];
				//cout << _i << ':' << P[_i] <<'?' << dP[_i] << 'l' << P[_i] - GetPi(*it, sV, sCita) <<endl;

				if (maxPe < fabs(dP[_i]))
				{
					maxPe = fabs(dP[_i]);
				}
			}
			
			if (maxPe < e)
			{
				KP = 0;

				if (KQ == 0)
				{
					break;
				}
			}else{
				B1.SolveEquation(dP);
				
				for (vector<double>::iterator it = cita.begin(); it != cita.end(); it++)
				{
					int _i = distance(cita.begin(), it);
					*it += dP[_i];
					sCita[B1Map[_i]] = *it;
					//cout << "cita" << B1Map[_i] << ':' << dP[_i] << endl ;
				}

				KQ = 1;
			}

			double maxQe = -1;
			for (vector<int>::const_iterator it = B2Map.begin(); it != B2Map.end(); it++)
			{	
				int _i = distance<vector<int>::const_iterator>(B2Map.begin(), it);
				dQ[_i] = (Q[_i] - GetQi(*it, sV, sCita)) / V[_i]; 
				//cout << Q[_i] - GetQi(*it, sV, sCita) << endl;

				if (maxQe < fabs(dQ[_i]))
				{
					maxQe = fabs(dQ[_i]);
				}
			}
			
			if (maxQe < e)
			{
				KQ = 0;

				if (KP == 0)
				{
					break;
				}
			}else{
				B2.SolveEquation(dQ);
				
				for (vector<double>::iterator it = V.begin(); it != V.end(); it++)
				{
					int _i = distance(V.begin(), it);
					*it += dQ[_i];
					//cout << "dV" << _i << ':' << dQ[_i]<< endl ;
				}

				KP = 1;
			}		
			k++;
			
		}
		Times = k;

		//Update the Bus[]
		for (vector<Pair>::const_iterator it = NodeIdMap.begin(); it != NodeIdMap.end(); it++)
		{
			int i = distance<vector<Pair>::const_iterator>(NodeIdMap.begin(), it);

			BusAndBranch::Bus *pBu = BAB.GetBus(it->j);

			pBu->V = sV[i];
			pBu->cita = sCita[i];
			pBu->P = GetPi(i, sV, sCita);
			pBu->Q = GetQi(i, sV, sCita);
		}

		BranchFlowCal();
	}

	void FastPowerFlow::Print() const
	{
		cout << "Bus Information:" << endl;
		BAB.PrintBuses();

		cout << "Branch apparent power:" << endl;

		for (int i = 0; i<(int)Pij.size(); i++)
		{
			cout << "S[" << Pij[i].i << "][" << Pij[i].j << "] = " << Pij[i].data << "+j" << Qij[i].data << endl;
		}
	}

	int FastPowerFlow::GetB2Virtual(int i) const
	{
		for (vector<int>::const_iterator it = B2Map.begin(); it != B2Map.end(); it++)
		{
			if (*it == i)
			{
				return distance(B2Map.begin(), it);
			}
		}
		return -1;

	}

	int FastPowerFlow::GetB1Virtual(int i) const
	{
		for (vector<int>::const_iterator it = B1Map.begin(); it != B1Map.end(); it++)
		{
			if (*it == i)
			{
				return distance(B1Map.begin(), it);
			}
		}
		return -1;
	}

	void FastPowerFlow::ToFile(string file)
	{
		fstream fp;
		fp.open(file.c_str(), ios::out | ios::app);

		string func;
		switch (type)
		{
		case XB:
			func = "XB:";
			break;
		case BX:
			func = "BX:";
			break;
		case PQ:
			func = "PQ:";
			break;
		default:
			throw _Exception("No this kind of algorithm(In ToFile)");
			break;
		}

		fp << func << endl;

		fp << "iterative times:" << Times << endl;

		fp << "Bus information:" << endl;
		for (int i=0; i<BAB.GetBusSize(); i++)
		{
			BusAndBranch::Bus const*  pBu = BAB.GetBus(i);
			fp << pBu->i << ":V=" << pBu->V << " cita=" << pBu->cita << "(" << pBu->cita / PI * 180<< ") S=" << pBu->P << "+j" << pBu->Q << endl;
		}

		fp << "Branch apparent power:" << endl;

		for (int i = 0; i<(int)Pij.size(); i++)
		{
			fp << "S[" << Pij[i].i << "][" << Pij[i].j << "] = " << Pij[i].data << "+j" << Qij[i].data << endl;
		}

		fp << "////////////////Copyright (c) Peng Yao//////////////\r\n\r\n";
		fp.flush();
		fp.clear();
		fp.close();
	}

	NRPowerFlow::NRPowerFlow()
	{
	}

	NRPowerFlow::~NRPowerFlow()
	{
	}

	void NRPowerFlow::CreateYMatrix()
	{
		SortBABByType();
		for (vector<Pair>::const_iterator it = NodeIdMap.begin(); it != NodeIdMap.end(); it++)
		{
			int realId = BAB.GetBus(it->j)->i;
			Y.push_back(ThreeElements<complex<double>>(it->i, it->i, complex<double>(0, BAB.GetBus(it->j)->b)));
			int major = Y.size() - 1;

			for (int i = 0; i < BAB.GetBranchSize(); i++)
			{
				BusAndBranch::Branch* pBr = BAB.GetBranch(i);

				if (pBr->i == realId || pBr->j == realId)
				{
					bool left = pBr->i == realId;
					complex<double> y = complex<double>(-1.0) / complex<double>(pBr->R , pBr->X), ym = complex<double>(-1.0) * y;
					switch (pBr->Type)
					{
					case BusAndBranch::Branch::LINE:
						break;
					case BusAndBranch::Branch::TRANSFORMER:
						y /= complex<double>(pBr->K);
						if (!left)
						{
							ym /= complex<double>(pow(pBr->K, 2.0));
						}
						break;
					default:
						throw _Exception("Unknown branch type(In CreateYMatrix)");
						break;
					}
					if (GetVirtualId(pBr->i) >= it->i && GetVirtualId(pBr->j) >= it->i)
					{
						Y.push_back(ThreeElements<complex<double>>(GetVirtualId(pBr->i), GetVirtualId(pBr->j), y));
						Y.push_back(ThreeElements<complex<double>>(GetVirtualId(pBr->j), GetVirtualId(pBr->i), y));
					}
					
					Y[major].data += ym + complex<double>(0.0, pBr->B);
				}
			}
		}

		int add = 0;
		if (nPQ + nPV > 2)
		{
			vector<ThreeElements<BYTE>> Yt;
			Yt.reserve(Y.size());
			int offset = nSlack;

			for (vector<ThreeElements<complex<double>>>::const_iterator it = Y.begin(); it != Y.end(); it++)
			{
				int i = it->i - offset, j = it->j - offset;

				if (i >= 0 && j >= 0)
				{
					Yt.push_back(ThreeElements<BYTE>(i, j, 0));
				}
			}
			Graph<BYTE> gh;
			gh.FromTEs(Yt);
			vector<int> ps;		
			vector<PowerSystem::Pair> psAdd;	
			gh.Optimization(ps, add, &psAdd);
			add *= 2;
			
		/*	gh.Print();
			PowerSystem::test_Print(ps);
			cout << "\n注入元个数" << add * 2<< endl;
			test_Print(psAdd);*/
			vector<Pair> temp;
			temp.push_back(NodeIdMap[0]);
			for (vector<int>::iterator it = ps.begin(); it != ps.end(); it++)
			{
				temp.push_back(NodeIdMap[(*it) + 1]);
			}
			NodeIdMap.assign(temp.begin(), temp.end());

			for (vector<Pair>::iterator it = NodeIdMap.begin(); it != NodeIdMap.end(); it++)
			{
				it->i = distance(NodeIdMap.begin(), it);
			}

			for (vector<ThreeElements<complex<double>>>::iterator it = Y.begin(); it != Y.end(); it++)
			{
				bool s1 = false, s2 = false;

				for (int i = 0; i<(int)ps.size(); i++)
				{
					if (!s1 && it->i == ps[i] + 1)
					{
						it->i = i + 1;
						s1 = true;
					}

					if (!s2 && it->j == ps[i] + 1)
					{
						it->j = i + 1;
						s2 = true;
					}

					if (s1 && s2)
					{
						break;
					}
				}
			}	
		}
		JacobiSize = 4 * (Y.size() + add);
		SparseMatrix<complex<double>>::SortTEs(Y, GetNodeNum());
		SparseMatrix<complex<double>> sm(Y, Y.size(), GetNodeNum());
		Matrix<complex<double>> m;
		sm.ToMatrix(m);
		m.ToFile("NY.txt");	
	}

	void NRPowerFlow::InitJacobi()
	{
		JacobiPos.reserve(GetNodeNum());
		JacobiDis.reserve(nPV + nPQ);
		int i = 0;
		for (vector<Pair>::const_iterator it = NodeIdMap.begin(); it != NodeIdMap.end(); it++)
		{
			BusAndBranch::Bus *pBu = BAB.GetBus(it->j);	

			switch (pBu->Type)
			{
			case BusAndBranch::Bus::PQ:
				JacobiPos.push_back(Pair(i, i + 1));
				i += 2;
				JacobiDis.push_back(2);
				break;
			case BusAndBranch::Bus::PV:
				JacobiPos.push_back(Pair(i++, -1));
				JacobiDis.push_back(1);
				break;
			case BusAndBranch::Bus::SLACK:
				JacobiPos.push_back(Pair(-1, -1));
				break;
			default:
				throw _Exception("Unknown branch type(In InitJacobiPostion)");
				break;
			}
		}

	}

	int NRPowerFlow::GetDistance(int i, int j)
	{
		if (i >= (int)JacobiDis.size() || j >= (int)JacobiDis.size() + 1 || i < 0 || j < 0)
		{
			throw _Exception("Index out of range of JacobiDis(In GetDistance)");
		}

		int ret = 0;
		for (int _i = i; _i < j; _i++)
		{
			ret += JacobiDis[_i];
		}

		return ret;
	}

	void NRPowerFlow::FlowCal(int type, int times, double e)
	{
		Pij.clear();
		Qij.clear();
		InitJacobi();
		int CalNodeNum = nPV + 2 * nPQ;
		vector<double> dS(CalNodeNum), Ssp(CalNodeNum), Si(CalNodeNum), dU(CalNodeNum), U(CalNodeNum);
		vector<double> sV(GetNodeNum()), sCita(GetNodeNum(), 0.0);

		int i = 0;
		for (vector<Pair>::const_iterator it = NodeIdMap.begin(); it != NodeIdMap.end(); it++)
		{	
			if (it->i - 1 >= 0)
			{
				BusAndBranch::Bus *pBu = BAB.GetBus(it->j);	
				
				switch (pBu->Type)
				{
				case BusAndBranch::Bus::PQ:	
					Ssp[i] = pBu->GP - pBu->LP;
					U[i++] = 0;
					Ssp[i] = pBu->GQ - pBu->LQ;	
					U[i++] = 1.0;
					break;
				case BusAndBranch::Bus::PV:
					Ssp[i] = pBu->GP - pBu->LP;
					U[i++] = 0.0;
					break;
				default:
					break;
				}
			}
		}

		int k = 0;
		while (k <= times)
		{
			for (vector<Pair>::const_iterator it = NodeIdMap.begin(); it != NodeIdMap.end(); it++)
			{
				BusAndBranch::Bus *pBu = BAB.GetBus(it->j);	

				switch (pBu->Type)
				{
				case BusAndBranch::Bus::PQ:
					sV[it->i] = U[JacobiPos[it->i].j];
					sCita[it->i] = U[JacobiPos[it->i].i];
					break;
				case BusAndBranch::Bus::PV:
					sV[it->i] = pBu->V;
					sCita[it->i] = U[JacobiPos[it->i].i];
					break;
				case BusAndBranch::Bus::SLACK:
					sV[it->i] = pBu->V;
					sCita[it->i] = pBu->cita;
					break;
				default:
					break;
				}
			}
			
			int i = 0;
			double max_e = -1;
			for (vector<Pair>::const_iterator it = NodeIdMap.begin(); it != NodeIdMap.end(); it++)
			{
				if (it->i - 1 >= 0)
				{
					BusAndBranch::Bus *pBu = BAB.GetBus(it->j);	
					
					switch (pBu->Type)
					{
					case BusAndBranch::Bus::PQ:	
						Si[i] = GetPi(it->i, sV, sCita);
						dS[i] = Ssp[i] - Si[i];
						if (max_e < fabs(dS[i]))
						{
							max_e = fabs(dS[i]);
						}
						i++;

						Si[i] = GetQi(it->i, sV, sCita);
						dS[i] = Ssp[i] - Si[i];
						if (max_e < fabs(dS[i]))
						{
							max_e = fabs(dS[i]);
						}
						i++;
						break;
					case BusAndBranch::Bus::PV:
						Si[i] = GetPi(it->i, sV, sCita);
						dS[i] = Ssp[i] - Si[i];
						if (max_e < fabs(dS[i]))
						{
							max_e = fabs(dS[i]);
						}
						i++;
						break;
					default:
						break;
					}
				}
			}

			//cout <<max_e << endl;
			if (max_e < e)
			{
				break;
			}

			vector<ThreeElements<double>> JacobiMatrix;
			JacobiMatrix.reserve(JacobiSize);
			for (vector<ThreeElements<complex<double>>>::const_iterator it = Y.begin(); it != Y.end(); it++)
			{
				int i = it->i - nSlack, j = it->j - nSlack;
				
				if (i >=0 && j >= 0)
				{
					BusAndBranch::Bus *pBuI = BAB.GetBus(GetSystemMap(it->i));
					if (i == j)
					{
						int offset = GetDistance(0, i);
						int _i = offset;
						double Citai = sCita[it->i];
						double Pi = Si[offset];					
						double Vi = sV[it->i];
						if (pBuI->Type == BusAndBranch::Bus::PQ)
						{	
							double Qi = Si[offset + 1];
							JacobiMatrix.push_back(ThreeElements<double>(_i, _i, pow(Vi, 2.0) * it->data.imag() + Qi));
							JacobiMatrix.push_back(ThreeElements<double>(_i, _i + 1, -pow(Vi, 2.0) * it->data.real() - Pi));
							JacobiMatrix.push_back(ThreeElements<double>(_i + 1, _i, pow(Vi, 2.0) * it->data.real() - Pi));
							JacobiMatrix.push_back(ThreeElements<double>(_i + 1, _i + 1, pow(Vi, 2.0) * it->data.imag() - Qi));
						}else{
							/* PV */
							JacobiMatrix.push_back(ThreeElements<double>(_i, _i, pow(Vi, 2.0) * it->data.imag() + GetQi(it->i, sV, sCita)));
						}
					}else{
						BusAndBranch::Bus *pBuJ = BAB.GetBus(GetSystemMap(it->j));
						int offset_i = GetDistance(0, i), offset_j = GetDistance(0, j);
						int _i = offset_i, _j = offset_j;
						double Citai = sCita[it->i];
						double Vi = sV[it->i];
						double Citaj = sCita[it->j];
						double Vj = sV[it->j];
						if (pBuI->Type == BusAndBranch::Bus::PQ)
						{
							double Hij = Vi * Vj * (it->data.imag() * cos(Citai - Citaj) - it->data.real() * sin(Citai - Citaj));
							double Mij = Vi * Vj * (it->data.real() * cos(Citai - Citaj) + it->data.imag() * sin(Citai - Citaj));
							JacobiMatrix.push_back(ThreeElements<double>(_i, _j, Hij));
							JacobiMatrix.push_back(ThreeElements<double>(_i + 1, _j, Mij));
							if (pBuJ->Type == BusAndBranch::Bus::PQ)
							{
								double Nij = -Mij, Lij = Hij;
								JacobiMatrix.push_back(ThreeElements<double>(_i, _j + 1, Nij));
								JacobiMatrix.push_back(ThreeElements<double>(_i + 1, _j + 1, Lij));
							}
						}else{
							/* PV */
							double Hij = Vi * Vj * (it->data.imag() * cos(Citai - Citaj) - it->data.real() * sin(Citai - Citaj));
							JacobiMatrix.push_back(ThreeElements<double>(_i, _j, Hij));
							if (pBuJ->Type == BusAndBranch::Bus::PQ)
							{
								double Nij = Vi * Vj * (it->data.real() * cos(Citai - Citaj) + it->data.imag() * sin(Citai - Citaj));
								JacobiMatrix.push_back(ThreeElements<double>(_i, _j + 1, Nij));
							}						
						}
					}
				}
			}

			int rc = GetDistance(0, nPQ + nPV);
			SparseMatrix<double>::SortTEs(JacobiMatrix, rc);
			SparseMatrix<double> sJacobi(JacobiMatrix, JacobiSize, rc);
			sJacobi.LDUDecomposition(false);
			sJacobi.SolveEquation(dS);
			i = 0;
			for (vector<Pair>::const_iterator it = NodeIdMap.begin(); it != NodeIdMap.end(); it++)
			{	
				if (it->i - 1 >= 0)
				{
					BusAndBranch::Bus *pBu = BAB.GetBus(it->j);	
					
					switch (pBu->Type)
					{
					case BusAndBranch::Bus::PQ:	
						U[i] -= dS[i];
						i++;
						U[i] -= dS[i] * U[i];
						i++;
						break;
					case BusAndBranch::Bus::PV:
						U[i] -= dS[i];
						i++;
						break;
					default:
						break;
					}
				}
			}
			k++;
		}
		Times = k;

		//Update the Bus[]
		for (vector<Pair>::const_iterator it = NodeIdMap.begin(); it != NodeIdMap.end(); it++)
		{
			int i = distance<vector<Pair>::const_iterator>(NodeIdMap.begin(), it);

			BusAndBranch::Bus *pBu = BAB.GetBus(it->j);

			pBu->V = sV[i];
			pBu->cita = sCita[i];
			pBu->P = GetPi(i, sV, sCita);
			pBu->Q = GetQi(i, sV, sCita);
		}

		BranchFlowCal();

	}

	void NRPowerFlow::Print() const
	{
		cout << "Bus Information:" << endl;
		BAB.PrintBuses();

		cout << "Branch apparent power:" << endl;

		for (int i = 0; i<(int)Pij.size(); i++)
		{
			cout << "S[" << Pij[i].i << "][" << Pij[i].j << "] = " << Pij[i].data << "+j" << Qij[i].data << endl;
		}
	}

	void NRPowerFlow::ToFile(string file)
	{
		fstream fp;
		fp.open(file.c_str(), ios::out | ios::app);

		string func = "NR:";
		fp << func << endl;

		fp << "Iterative times:" << Times << endl;

		fp << "Bus information:" << endl;
		for (int i=0; i<BAB.GetBusSize(); i++)
		{
			BusAndBranch::Bus const*  pBu = BAB.GetBus(i);
			fp << pBu->i << ":V=" << pBu->V << " cita=" << pBu->cita << "(" << pBu->cita / PI * 180<< ") S=" << pBu->P << "+j" << pBu->Q << endl;
		}

		fp << "Branch apparent power:" << endl;

		for (int i = 0; i<(int)Pij.size(); i++)
		{
			fp << "S[" << Pij[i].i << "][" << Pij[i].j << "] = " << Pij[i].data << "+j" << Qij[i].data << endl;
		}
		fp << "////////////////Copyright (c) Peng Yao//////////////\r\n\r\n";
		fp.flush();
		fp.clear();
		fp.close();
	}
}


