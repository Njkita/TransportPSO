#include <iostream>
#include <random>
#include <math.h>
#include <iomanip>
#include <utility>
#include <algorithm>
#include "RoadMap.h"
#include "AstarPSOfinal.h"

double Dijkstra(Graph GR, int st, int en)
{
	int index, u, m = st + 1;
	bool visited[quantity];
	double distance[quantity];

	for (int i = 0; i < quantity; i++)
	{
		distance[i] = INT_MAX; 
		visited[i] = false;
	}

	distance[st] = 0;

	for (int count = 0; count < quantity - 1; count++)
	{
		int min = INT_MAX;

		for (int i = 0; i < quantity; i++)
			if (!visited[i] && distance[i] <= min)
			{
				min = distance[i]; 
				index = i;
			}

		u = index;
		visited[u] = true;

		for (int i = 0; i < quantity; i++)
			if (!visited[i] && GR.roads[u][i] && distance[u] != INT_MAX && distance[u] + GR.roads[u][i] < distance[i])
				distance[i] = distance[u] + GR.roads[u][i];
	}

	

	return distance[en];
}

//Алгоритм PSO
//--------------------------------------Параметры агентов\шахт\терминалов---------------------------------------------------

// Количество агентов (например, поездов)
const int N = 6;
// Количество задач (шахт)
const int M = mine;
// Размер популяции частиц
const int P = 6 * M;
// Коэффициент инерции
const double omega = 0.7;
// Вес локального лучшего решения
const double C1 = 1;
// Вес глобального лучшего решения
const double C2 = 1;

//Оптимизирующая функция
double fitness(Graph GR, int tasks[M], int agents[M])
{

	double g1 = 1, g2 = 1, g3 = 1;

	//запись информации о терминалах в файл
	std::ofstream out;
	out.open("terminal.txt");

	out << quantity << " " << mine << " " << terminal << " " << N << std::endl;

	for (int i = 0; i < M; i++)
		for (int j = 0; j < M; j++)
			if (agents[j] == i)
				out << agents[j] << " ";

	out << std::endl;

	for (int i = 0; i < M; i++)
		for (int j = 0; j < M; j++)
			if (agents[j] == i)
				out << j << " ";

	out << std::endl;

	for (int i = 0; i < M; i++)
		for (int j = 0; j < M; j++)
			if (agents[j] == i)
				out << tasks[j] + mine << " ";

	out << std::endl;

	out.close();

	//Массив для хранения времени работы агентов
	std::pair<int, double> times[N]{};
	for (int i = 0; i < N; i++) { times[i].first = -1; }

	//Массивы для отметки посещенных шахт и терминалов
	bool Mines[mine]{};
	bool Terminals[terminal]{};

	//Время работы агентов
	//i - номер агента
	//tasks[i] - номер терминала, соответствующий агенту
	double z1 = 0;

	//Доля невывезенного груза из шахт
	double z2 = 0;

	//Доля недовезенного груза до терминалов
	double z3 = 0;

	//Массив, содержащий количество груза, которое агенты везут на терминал
	double cargo[N]{};

	//Временные массивы для хранения потребностей шахт и терминалов
	double needMines[mine]{};
	double needTerminals[terminal]{};

	for (int i = 0; i < mine; i++)
	{
		needMines[i] = GR.wantMine[i];
	}

	for (int i = 0; i < terminal; i++)
	{
		needTerminals[i] = GR.wantTerminal[i];
	}

	//Запуск симуляции маршрутов
	run_simulation("Graph3.txt", "terminal.txt");

	//Чтение времени работы агентов из файла
	std::ifstream in;
	in.open("pso.txt");

	for (int i = 0; i < N; i++)
		in >> times[i].second;

	in.close();

	for (int i = 0; i < M; i++)
	{
		Mines[i] = true;
		Terminals[tasks[i]] = true;

		/*
		if (times[agents[i]].first == -1)
			times[agents[i]].second = Dijkstra(GR, i, tasks[i]);
		else
			times[agents[i]].second += Dijkstra(GR, times[agents[i]].first, i) + Dijkstra(GR, i, tasks[i]);
			*/

		//Количество вывезенного груза
		cargo[agents[i]] += needMines[i];

		//Обновление остатка груза
		needMines[i] = 0;

		//Разгрузка груза на терминале
		if (cargo[agents[i]] >= needTerminals[tasks[i]])
		{
			cargo[agents[i]] -= needTerminals[tasks[i]];
			needTerminals[tasks[i]] = 0;
		}
		else
		{
			needTerminals[tasks[i]] -= cargo[agents[i]];
			cargo[agents[i]] = 0;
		}

		times[agents[i]].first = tasks[i];
	}

	//Доля недовезенного груза до терминалов
	for (int i = 0; i < terminal; i++)
	{
		if (needTerminals[i] > 0)
			z3 += needTerminals[i] / GR.wantTerminal[i];
	}

	//Доля невывезенного груза из шахт
	for (int i = 0; i < mine; i++)
	{
		if (needMines[i] > 0)
			z2 += needMines[i] / GR.wantMine[i];
    }

	//Время работы, чем больше агентов простаивает, тем хуже
	for (int i = 0; i < N; i++)
		if (times[i].first == -1)
			z1 += N;
		else
			z1 += times[i].second / 24;

	//Доля груза, который агенты везут, но не разгрузили
	for (int i = 0; i < N; i++)
		if (cargo[i] > 0)
			z2 += cargo[i];

	//Штраф за шахты и терминалы, которые не были посещены
	for (int i = 0; i < mine; i++)
	{
		if (Mines[i] == false)
			z2 += mine;
	}

	for (int i = 0; i < terminal; i++)
	{
		if (Terminals[i] == false)
			z3 += terminal;
	}

	//Штраф за неравномерное распределение задач между агентами (чтобы избежать ситуации, когда один агент выполняет все задачи, а другие простаивают)
	std::vector<int> counts(N, 0);

	for (int i = 0; i < M; i++)
		counts[agents[i]]++;

	if (*std::max_element(begin(counts), end(counts)) > std::ceil((double)M / N))
		z1 += N;

	//Нормализация значений на диапазоны
	z1 /= N;
	z2 /= mine;
	z3 /= terminal;

	return g1 * z1 + g2 * z2 + g3 * z3;
}

//Функция для вывода текущего состояния грузов, на основе данных из файла
double inspection(Graph GR, int tasks[M], int agents[M])
{
	//Массив, содержащий количество груза, которое агенты везут на терминал
	double cargo[N]{};

	for (int i = 0; i < M; i++)
	{
		//Количество вывезенного груза
		cargo[agents[i]] += GR.wantMine[i];

		//Обновление остатка груза
		GR.wantMine[i] = 0;

		//Разгрузка груза на терминале
		if (cargo[agents[i]] >= GR.wantTerminal[tasks[i]])
		{
			cargo[agents[i]] -= GR.wantTerminal[tasks[i]];
			GR.wantTerminal[tasks[i]] = 0;
		}
		else
		{
			GR.wantTerminal[tasks[i]] -= cargo[agents[i]];
			cargo[agents[i]] = 0;
		}
	}

	std::cout << std::endl;

	for (int i = 0; i < N; i++)
		std::cout << i << "\t";

	std::cout << std::endl;

	for (int i = 0; i < N; i++)
		std::cout << cargo[i] << "\t";

	std::cout << std::endl;

	return 0;
}

double r()
{
	std::random_device device;

	std::mt19937_64 engine(device());

	std::uniform_real_distribution<> distribution(0.0, 1.0);

	const auto U = distribution(engine);

	return U;
}

int initPSO(int max)
{
	return floor(max * r());
}

double initV(int max)
{
	double sign = r();

	if (sign >= 0.5)
		return floor(max * r());
	else
		return -floor(max * r());
}

int main()
{
	Graph GR("Graph.txt");
	GR.print();
	std::cout << Dijkstra(GR, 0, 10) << std::endl;

	//Массив, содержащий задачи и агентов, в первом измерении хранятся частицы, во втором - задачи
	//Каждая частица представляет собой возможное решение
	int Tasks[P][M]{};
	int Agents[P][M]{};

	//Скорости частиц
	double vTasks[P][M]{};
	double vAgents[P][M]{};

	//Индекс лучшей глобальной частицы
	int Gbest = 0;

	//Массив лучших позиций частиц
	int PBestTasks[P][M]{};
	int PBestAgents[P][M]{};

	//Массив значений функции приспособленности для лучших позиций частиц
	double PBestFitness[P];

	//Этап 1: Инициализация
	//Инициализация позиций
	for (int i = 0; i < P; i++)
		for (int j = 0; j < M; j++)
		{
			Tasks[i][j]  = initPSO(terminal);
			Agents[i][j] = initPSO(N);
		}
	
	//Инициализация скоростей
	for (int i = 0; i < P; i++)
		for (int j = 0; j < M; j++)
		{
			vTasks[i][j]  = initV(terminal);
			vAgents[i][j] = initV(N);
		}

	//Инициализация лучших позиций частиц
	for (int i = 0; i < P; i++)
		for (int j = 0; j < M; j++)
		{
			PBestTasks[i][j]  = Tasks[i][j];
			PBestAgents[i][j] = Agents[i][j];
		}

#if 1
	//Вывод информации о начальной популяции:
	std::cout << std::setprecision(1);
	std::cout << std::fixed;

	std::cout << "--------------------------" << std::endl;
	std::cout << "For generation [" << 0 << "]" << std::endl;
	std::cout << "--------------------------" << std::endl;

	for (int i = 0; i < P; i++)
	{
		for (int j = 0; j < M; j++)
			std::cout << std::setw(6) << vTasks[i][j];

		std::cout << "\t";

		for (int j = 0; j < M; j++)
			std::cout << std::setw(5) << Tasks[i][j];

		std::cout << "\t";

		std::cout << std::setw(5) << fitness(GR, Tasks[i], Agents[i]) << std::endl;

		for (int j = 0; j < M; j++)
			std::cout << std::setw(6) << vAgents[i][j];

		std::cout << "\t";

		for (int j = 0; j < M; j++)
			std::cout << std::setw(5) << Agents[i][j];

		std::cout << std::endl << std::endl;
	}
#endif
	/*
	//Вывод значений функции приспособленности для начальной популяции
	std::cout << std::setprecision(1);
	std::cout << std::fixed;

	for (int i = 0; i < P; i++)
		std::cout << std::setw(10) << fitness(GR, Tasks[i], Agents[i]) << "\t";
	
	std::cout << std::endl << std::endl;
	*/
	//Этап 2: Оптимизация с помощью PSO
	for (int l = 0; l < 36; l++)
	{
		//Поиск Gbest
		double min = fitness(GR, Tasks[0], Agents[0]);
		Gbest = 0;

		for (int i = 1; i < P; i++)
		{
			if (min > fitness(GR, Tasks[i], Agents[i]))
			{
				min = fitness(GR, Tasks[i], Agents[i]);
				Gbest = i;
			}
		}

		//Обновление лучших позиций частиц
		for (int i = 0; i < P; i++)
		{
			if (fitness(GR, Tasks[i], Agents[i]) < fitness(GR, PBestTasks[i], PBestAgents[i]))
				for (int j = 0; j < M; j++)
				{
					PBestTasks[i][j] = Tasks[i][j];
					PBestAgents[i][j] = Agents[i][j];
				}
		}

		//Обновление скоростей и позиций частиц
		//Скорости:
		for (int i = 0; i < P; i++)
			for (int j = 0; j < M; j++)
			{
				vTasks[i][j] = omega * vTasks[i][j] + C1 * r() * (PBestTasks[i][j] - Tasks[i][j]) + C2 * r() * (Tasks[Gbest][j] - Tasks[i][j]);
				vAgents[i][j] = omega * vAgents[i][j] + C1 * r() * (PBestAgents[i][j] - Agents[i][j]) + C2 * r() * (Agents[Gbest][j] - Agents[i][j]);
			}

		//Позиции:
		for (int i = 0; i < P; i++)
			for (int j = 0; j < M; j++)
			{
				Tasks[i][j] = floor(abs(Tasks[i][j] + vTasks[i][j]));

				if (Tasks[i][j] >= terminal)
					Tasks[i][j] = terminal - 1;

				Agents[i][j] = floor(abs(Agents[i][j] + vAgents[i][j]));

				if (Agents[i][j] >= N)
					Agents[i][j] = N - 1;
			}

#if 1
		//Вывод информации о текущей популяции:
		std::cout << "--------------------------" << std::endl;
		std::cout << "For generation [" << l + 1 << "]" << std::endl;
		std::cout << "--------------------------" << std::endl;

		for (int i = 0; i < P; i++)
		{
			for (int j = 0; j < M; j++)
				std::cout << std::setw(6) << vTasks[i][j];

			std::cout << "\t";

			for (int j = 0; j < M; j++)
				std::cout << std::setw(5) << Tasks[i][j];

			std::cout << "\t";

			std::cout << std::setw(5) << fitness(GR, Tasks[i], Agents[i]) << std::endl;

			for (int j = 0; j < M; j++)
				std::cout << std::setw(6) << vAgents[i][j];

			std::cout << "\t";

			for (int j = 0; j < M; j++)
				std::cout << std::setw(5) << Agents[i][j];

			std::cout << std::endl << std::endl;
		}
#endif
		/*
		//Вывод значений функции приспособленности для текущей популяции
		for (int i = 0; i < P; i++)
			std::cout << std::setw(10) << fitness(GR, Tasks[i], Agents[i]) << "\t";

		std::cout << std::endl << std::endl;
		*/
	}

	//Финальное обновление лучших позиций частиц
	for (int i = 0; i < P; i++)
	{
		if (fitness(GR, Tasks[i], Agents[i]) < fitness(GR, PBestTasks[i], PBestAgents[i]))
			for (int j = 0; j < M; j++)
			{
				PBestTasks[i][j] = Tasks[i][j];
				PBestAgents[i][j] = Agents[i][j];
			}
	}

	//Поиск лучшего решения среди всех частиц:
	double min = fitness(GR, PBestTasks[0], PBestAgents[0]);
	Gbest = 0;

	for (int i = 1; i < P; i++)
	{
		if (min > fitness(GR, PBestTasks[i], PBestAgents[i]))
		{
			min = fitness(GR, PBestTasks[i], PBestAgents[i]);
			Gbest = i;
		}
	}

//	std::cout << "0 0 1 2 3 3 4 5 5" << std::endl;
	for (int j = 0; j < M; j++)
		std::cout << PBestTasks[Gbest][j] << " ";

	std::cout << std::endl;

	//Вывод терминалов, которые не были задействованы в лучшем решении (если есть такие термины, нужно проверить решение)
	for (int i = 0; i < terminal; i++)
	{
		bool flag = false;

		for (int j = 0; j < M; j++)
			if (PBestTasks[Gbest][j] == i)
				flag = true;

		if (flag == false)
			std::cout << i << "\t";
	}

	std::cout << std::endl;

	for (int j = 0; j < M; j++)
		std::cout << PBestAgents[Gbest][j] << " ";

	std::cout << std::endl;

	//Вывод агентов, которые не были задействованы в лучшем решении (если есть такие агенты, нужно проверить решение)
	for (int i = 0; i < N; i++)
	{
		bool flag = false;

		for (int j = 0; j < M; j++)
			if (PBestAgents[Gbest][j] == i)
				flag = true;

		if (flag == false)
			std::cout << i << "\t";
	}

	inspection(GR, PBestTasks[Gbest], PBestAgents[Gbest]);

	std::cout << "::::::::::::::::::::::::::" << std::endl;
	std::cout << fitness(GR, PBestTasks[Gbest], PBestAgents[Gbest]) << std::endl;
}