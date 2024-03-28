#include <iostream>
#include <Windows.h>
#include<iomanip>

const int Nsize = 2000;// 10000;//120 500 1400 for 64KB(128KB) 2*256KB 3MB
int vec[Nsize]{};
int mat[Nsize][Nsize]{};
int sum[Nsize]{};

void Initial(int Nsize) {
    for (int i = 0; i < Nsize; i++)
    {
        for (int j = i; j < Nsize; j++)
        {
            mat[i][j] = 1;
        }
    }
    for (int i = 0; i < Nsize; i++)
    {
        vec[i] = 1;
    }
    for (int j = 0; j < Nsize; j++)
    {
        sum[j] = 0;
    }
    // check �����Ƿ���ȷ��ʼ��
    //std::cout << "��ʼ����:" << std::endl;
    //for (int i = 0; i < Nsize; i++)
    //{
    //    for (int j = 0; j < Nsize; j++) {
    //        std::cout << mat[i][j];
    //        if (j < Nsize - 1)
    //            std::cout << ' ';
    //        else
    //            std::cout << std::endl;
    //    }
    //}
}

// ���������
void Zero_Sum(int Nsize) {
    for (int j = 0; j < Nsize; j++)
    {
        sum[j] = 0;
    }
}

// ���������
void Check_Sum(int Nsize, int repeat) {
    for (int i = 0; i < Nsize; i++)
    {
        std::cout << sum[i] / repeat;
        if (i < Nsize - 1)
            std::cout << ' ';
        else
            std::cout << std::endl;
    }
}

// ƽ���㷨���������
void Col_Calc(int Nsize) {
    for (int j = 0; j < Nsize; j++)
    {
        for (int i = 0; i < Nsize; i++)
        {
            sum[j] += vec[i] * mat[i][j];
        }
    }
}

// Cache�㷨���������
void Row_Calc(int Nsize) {
    for (int i = 0; i < Nsize; i++)
    {
        for (int j = 0; j < Nsize; j++)
        {
            sum[j] += vec[i] * mat[i][j];
        }
    }
}



// Test01 �������ݹ�ģ�����ܵ�Ӱ��
void Size_Test() {
    // �������ݹ�ģ�仯�Գ������е�Ӱ��
    // ���ַ�����ƽ������ʱ��

    int min_test_size = 10;// 100 -> 10000 // 10 -> 2000
    int max_test_size = 2000;
    int test_step = 10;// 50;// ���ı䲽��
    LARGE_INTEGER frequency, t1, t2;
    int repeat = 0;// ����������ʱ��϶̵ĳ������ѭ������
    double elapsedTime = 0;

    // ��ȡʱ��Ƶ��
    if (!QueryPerformanceFrequency(&frequency)) {
        std::cout << "QueryPerformanceFrequency failed!\n";
        return;
    }

    // ���в��Բ�������
    for (int n = min_test_size; n <= max_test_size; n+=test_step)
    {
        // ���в��� �м���
        Zero_Sum(max_test_size);
        QueryPerformanceCounter(&t1);
        for (repeat = 1; ; repeat++)// Col_Calc
        {
            Col_Calc(n);
            QueryPerformanceCounter(&t2);
            if ((t2.QuadPart - t1.QuadPart) * 1000.0 / frequency.QuadPart > 100)// ѭ��100ms
            {
                break;
            }
        }
        elapsedTime = (t2.QuadPart - t1.QuadPart) * 1000.0 / frequency.QuadPart / repeat; // ��λΪ����
        std::cout << "�м��� ��ģ: " << std::setw(4) << n << " ѭ������: " << std::setw(6) << repeat << " ƽ����ʱ: " << elapsedTime << " ms\n";

        // ���в��� �м���
        Zero_Sum(max_test_size);
        QueryPerformanceCounter(&t1);
        for (repeat = 1; ; repeat++)// Row_Calc
        {
            Row_Calc(n);
            QueryPerformanceCounter(&t2);
            if ((t2.QuadPart - t1.QuadPart) * 1000.0 / frequency.QuadPart > 100)// ѭ��100ms
            {
                break;
            }
        }
        elapsedTime = (t2.QuadPart - t1.QuadPart) * 1000.0 / frequency.QuadPart / repeat; // ��λΪ����
        std::cout << "�м��� ��ģ: " << std::setw(4) << n << " ѭ������: " << std::setw(6) << repeat << " ƽ����ʱ: " << elapsedTime << " ms\n";
    }

}

// Test02 ����CacheLine��С��ͨ���ֶ�����Row_Calc�Ĵ�С
void CacheLine_Size_Test() {

    int min_test_size = 2;// 2 -> 200
    int max_test_size = 200;
    int test_step = 2;// ����������
    LARGE_INTEGER frequency, t1, t2;
    int repeat = 0;// ����������ʱ��϶̵ĳ������ѭ������
    double elapsedTime = 0;

    // ��ȡʱ��Ƶ��
    if (!QueryPerformanceFrequency(&frequency)) {
        std::cout << "QueryPerformanceFrequency failed!\n";
        return;
    }

    // ���в��Բ�������
    for (int n = min_test_size; n <= max_test_size; n += test_step)
    {
        // ���в��� �м���
        Zero_Sum(max_test_size);
        QueryPerformanceCounter(&t1);
        for (repeat = 1; ; repeat++)// Col_Calc
        {
            Col_Calc(n);
            QueryPerformanceCounter(&t2);
            if ((t2.QuadPart - t1.QuadPart) * 1000.0 / frequency.QuadPart > 100)// ѭ��100ms
            {
                break;
            }
        }
        elapsedTime = (t2.QuadPart - t1.QuadPart) * 1000.0 / frequency.QuadPart / repeat; // ��λΪ����
        std::cout << "�м��� ��ģ: " << std::setw(4) << n << " ѭ������: " << std::setw(6) << repeat << " ƽ����ʱ: " << elapsedTime << " ms\n";

        // ���в��� �м���
        Zero_Sum(max_test_size);
        QueryPerformanceCounter(&t1);
        for (repeat = 1; ; repeat++)// Row_Calc
        {
            Row_Calc(n);
            QueryPerformanceCounter(&t2);
            if ((t2.QuadPart - t1.QuadPart) * 1000.0 / frequency.QuadPart > 100)// ѭ��100ms
            {
                break;
            }
        }
        elapsedTime = (t2.QuadPart - t1.QuadPart) * 1000.0 / frequency.QuadPart / repeat; // ��λΪ����
        std::cout << "�м��� ��ģ: " << std::setw(4) << n << " ѭ������: " << std::setw(6) << repeat << " ƽ����ʱ: " << elapsedTime << " ms\n";

        // ����������
    }

}

// Test03 ����Vtune����profileǰ 200-800 �� 800-1400 ��500��1100�������в��� �ظ�CacheLine���Է��� ����֤�����������½���L2��L3 Cache���������
void Pre_Vtune_Test() {

    int min_test_size = 800;// 800 -> 1400 // 200 -> 800
    int max_test_size = 1400;
    int test_step = 8;// ����������
    LARGE_INTEGER frequency, t1, t2;
    int repeat = 0;// ����������ʱ��϶̵ĳ������ѭ������
    double elapsedTime = 0;

    // ��ȡʱ��Ƶ��
    if (!QueryPerformanceFrequency(&frequency)) {
        std::cout << "QueryPerformanceFrequency failed!\n";
        return;
    }

    // ���в��Բ�������
    for (int n = min_test_size; n <= max_test_size; n += test_step)
    {
        // ���в��� �м���
        Zero_Sum(max_test_size);
        QueryPerformanceCounter(&t1);
        for (repeat = 1; ; repeat++)// Col_Calc
        {
            Col_Calc(n);
            QueryPerformanceCounter(&t2);
            if ((t2.QuadPart - t1.QuadPart) * 1000.0 / frequency.QuadPart > 100)// ѭ��100ms
            {
                break;
            }
        }
        elapsedTime = (t2.QuadPart - t1.QuadPart) * 1000.0 / frequency.QuadPart / repeat; // ��λΪ����
        std::cout << "Col_Calc " << std::setw(4) << n << " ѭ������: " << std::setw(6) << repeat << " ƽ����ʱ: " << elapsedTime << " ms\n";

        // ���в��� �м���
        Zero_Sum(max_test_size);
        QueryPerformanceCounter(&t1);
        for (repeat = 1; ; repeat++)// Row_Calc
        {
            Row_Calc(n);
            QueryPerformanceCounter(&t2);
            if ((t2.QuadPart - t1.QuadPart) * 1000.0 / frequency.QuadPart > 100)// ѭ��100ms
            {
                break;
            }
        }
        elapsedTime = (t2.QuadPart - t1.QuadPart) * 1000.0 / frequency.QuadPart / repeat; // ��λΪ����
        std::cout << "Row_Calc " << std::setw(4) << n << " ѭ������: " << std::setw(6) << repeat << " ƽ����ʱ: " << elapsedTime << " ms\n";

        // ����������
    }

}

// Test04 ��900 1100 1300���в��ԣ����Է���L3_Miss��Ϊ���� // רע��һ����һ��ģ
void Vtune_Test() {

    int min_test_size = 900;// 900 -> 1300
    int max_test_size = 1300;// ѡȡ���ʵ�������λ
    int test_step = 200;// ����������
    int n = min_test_size + test_step * 1;// �޸Ĺ�ģ ע�͵�����һ������

    LARGE_INTEGER frequency, t1, t2;
    int repeat = 0;// ����������ʱ��϶̵ĳ������ѭ������
    double elapsedTime = 0;

    // ��ȡʱ��Ƶ��
    if (!QueryPerformanceFrequency(&frequency)) {
        std::cout << "QueryPerformanceFrequency failed!\n";
        return;
    }

    // ���в��� �м���
    Zero_Sum(max_test_size);
    QueryPerformanceCounter(&t1);
    for (repeat = 1; ; repeat++)// Col_Calc
    {
        Col_Calc(n);
        QueryPerformanceCounter(&t2);
        if ((t2.QuadPart - t1.QuadPart) * 1000.0 / frequency.QuadPart > 100)// ѭ��100ms
        {
            break;
        }
    }
    elapsedTime = (t2.QuadPart - t1.QuadPart) * 1000.0 / frequency.QuadPart / repeat; // ��λΪ����
    std::cout << "�м��� ��ģ: " << std::setw(4) << n << " ѭ������: " << std::setw(6) << repeat << " ƽ����ʱ: " << elapsedTime << " ms\n";

    // ���в��� �м���
//    Zero_Sum(max_test_size);
//    QueryPerformanceCounter(&t1);
//    for (repeat = 1; ; repeat++)// Row_Calc
//    {
//        Row_Calc(n);
//        QueryPerformanceCounter(&t2);
//        if ((t2.QuadPart - t1.QuadPart) * 1000.0 / frequency.QuadPart > 100)// ѭ��100ms
//        {
//            break;
//        }
//    }
//    elapsedTime = (t2.QuadPart - t1.QuadPart) * 1000.0 / frequency.QuadPart / repeat; // ��λΪ����
//    std::cout << "�м��� ��ģ: " << std::setw(4) << n << " ѭ������: " << std::setw(6) << repeat << " ƽ����ʱ: " << elapsedTime << " ms\n";

}

int main()
{
    // Initial Nsize = 2000;
    Initial(Nsize);

    // �������� ��ģ�����ܵ�Ӱ�� ����Cache��С ����Cache Line����
    // Size_Test();

    // ̽��Cache Line��С(N = 2 -> 200)
    // CacheLine_Size_Test();

    // ����Vtune����profileǰ ��1100����(N = 800 -> 1400)���в��� �ظ�CacheLine���Է��� ����֤�����������½���L3 Cache���������
    // Pre_Vtune_Test();

    // ���Vtune����������Cache�������ϵĲ��� ֤�������ʼ�L3��С�����Ӱ���������ٶȣ�������ȡӰ�����ٶ�
    Vtune_Test();


    return 0;
}
