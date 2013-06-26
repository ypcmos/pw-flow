#include "pw_flow.h"

void main()
{
	try
	{	
		PowerSystem::FastPowerFlow fpf;
		fpf.InitBABFromFile("test14.txt");
		fpf.CreateYMatrix();
		
		fpf.FlowCal(PowerSystem::FastPowerFlow::XB , 100);
		fpf.ToFile();

		fpf.FlowCal(PowerSystem::FastPowerFlow::BX, 100);
		fpf.ToFile();

		fpf.FlowCal(PowerSystem::FastPowerFlow::PQ, 1000);
		fpf.ToFile();

		PowerSystem::NRPowerFlow nr;
		nr.InitBABFromFile("test14.txt");
		nr.CreateYMatrix();
		nr.FlowCal(1, 1000);
		nr.ToFile();
	}
	catch (const _Exception &e)
	{
		cout << e.Error() << endl;
	}
	catch(...)
	{
		cout << "UNKNOWN ERROR \r\n";
	}
	system("PAUSE");
}