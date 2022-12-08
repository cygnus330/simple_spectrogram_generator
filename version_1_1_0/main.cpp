#pragma GCC optimize("O3")
#pragma GCC optimize("Ofast")
#pragma GCC optimize("unroll-loops")
#pragma GCC target("avx,avx2,fma")

#include <bits/stdc++.h>
#include <Windows.h>
#include <chrono>
#include <cstdint>
#define endl "\n"

using namespace std;

struct WAV
{
    typedef unsigned char byte;
    byte ChunkID[4];      //"RIFF"
    byte ChunkSize[4];    //*
    byte Format[4];       //"WAVE"
    byte Subchunk1ID[4];  //"fmt "
    byte Subchunk1Size[4];//*
    byte AudioFormat[2];  //*, 1:pcm(pulse code modulation)
    byte NumChannels[2];  //*
    byte SampleRate[4];   //*
    byte ByteRate[4];     //*
    byte BlockAlign[2];   //*
    byte BitsPerSample[2];//*
    byte Subchunk2ID[4];  //"data"
    byte Subchunk2Size[4];//*
};
struct WAVt
{
    typedef unsigned char byte;
    byte ChunkID[4];    //"RIFF"
    int ChunkSize;      //*
    byte Format[4];     //"WAVE"
    byte Subchunk1ID[4];//"fmt "
    int Subchunk1Size;  //*
    short AudioFormat;  //*, 1:pcm(pulse code modulation)
    short NumChannels;  //*
    int SampleRate;     //*
    int ByteRate;       //*
    short BlockAlign;   //*
    short BitsPerSample;//*
    byte Subchunk2ID[4];//"data"
    int Subchunk2Size;  //*
};
int ler(unsigned char x[], int n)
{
    int t=0, s=4-n;
    for(int i=0; i<n; i++)
        t |= (x[i]<<((i+s)*8));
    return t>>(s*8);
}
WAVt readheader(WAV x)
{
    WAVt t;
    for(int i=0; i<4; i++) t.ChunkID[i]=x.ChunkID[i];
    t.ChunkSize=ler(x.ChunkSize, 4);
    for(int i=0; i<4; i++) t.Format[i]=x.Format[i];
    for(int i=0; i<4; i++) t.Subchunk1ID[i]=x.Subchunk1ID[i];
    t.Subchunk1Size=ler(x.Subchunk1Size, 4);
    t.AudioFormat=ler(x.AudioFormat, 2);
    t.NumChannels=ler(x.NumChannels, 2);
    t.SampleRate=ler(x.SampleRate, 4);
    t.ByteRate=ler(x.ByteRate, 4);
    t.BlockAlign=ler(x.BlockAlign, 2);
    t.BitsPerSample=ler(x.BitsPerSample, 2);
    for(int i=0; i<4; i++) t.Subchunk2ID[i]=x.Subchunk2ID[i];
    t.Subchunk2Size=ler(x.Subchunk2Size, 4);
    return t;
}
void lew(unsigned char *x, int t, int n)
{
    for(int i=0; i<n; i++, t>>=8)
        x[i] = t&0xFF;
}
WAV writeheader(WAVt x)
{
    WAV t;
    for(int i=0; i<4; i++) t.ChunkID[i]=x.ChunkID[i];
    lew(t.ChunkSize, x.ChunkSize, 4);
    for(int i=0; i<4; i++) t.Format[i]=x.Format[i];
    for(int i=0; i<4; i++) t.Subchunk1ID[i]=x.Subchunk1ID[i];
    lew(t.Subchunk1Size, x.Subchunk1Size, 4);
    lew(t.AudioFormat, x.AudioFormat, 2);
    lew(t.NumChannels, x.NumChannels, 2);
    lew(t.SampleRate, x.SampleRate, 4);
    lew(t.ByteRate, x.ByteRate, 4);
    lew(t.BlockAlign, x.BlockAlign, 2);
    lew(t.BitsPerSample, x.BitsPerSample, 2);
    for(int i=0; i<4; i++) t.Subchunk2ID[i]=x.Subchunk2ID[i];
    lew(t.Subchunk2Size, x.Subchunk2Size, 4);
    return t;
}

const double PI = 3.141592653589793;
const double E = 2.718281828459045;
int lengthAudio, channelsAudio, bytesAudio, samplerateAudio, minFreq = 0, maxFreq = 47;
string pitchTable[12] = {"C", "C#", "D", "D#", "E", "F", "F#", "G", "G#", "A", "A#", "B"};
unsigned char BUFFER;
vector<unsigned char> imageBUFFER(54);
vector<vector<vector<unsigned char>>> imageOriginal;

// text dot (3x5), abcdefg#
// 0A 1B 2C 3D 4E 5F 6G 7#
string txtDot[8][5] = {
    {" #  ", "# # ", "### ", "# # ", "# # "},
    {"##  ", "# # ", "##  ", "# # ", "##  "},
    {" ## ", "#   ", "#   ", "#   ", " ## "},
    {"##  ", "# # ", "# # ", "# # ", "##  "},
    {"### ", "#   ", "##  ", "#   ", "### "},
    {"### ", "#   ", "##  ", "#   ", "#   "},
    {" ## ", "#   ", "# # ", "# # ", " ## "},
    {"# # ", "### ", "# # ", "### ", "# # "}
};

// number dot (3x5)
string numDot[10][5] = {
    {"### ","# # ","# # ","# # ","### "},
    {" #  ","##  "," #  "," #  ","### "},
    {"### ","  # ","### ","#   ","### "},
    {"### ","  # "," ## ","  # ","### "},
    {"# # ","# # ","### ","  # ","  # "},
    {"### ","#   ","### ","  # ","### "},
    {"### ","#   ","### ","# # ","### "},
    {"### ","  # ","  # ","  # ","  # "},
    {"### ","# # ","### ","# # ","### "},
    {"### ","# # ","### ","  # ","### "}
};

// do DFT
double calcFourier(const vector<double> & PCM, int time, double frequency)
{
    int sampleCount = 512, cnt = 0;
    complex<double> calc(0, 0);
    for(int i = 0; i < sampleCount; i++)
    {
        if(time + i >= PCM.size() - 1)
            break;
        
        cnt++;
        // DFT in hwp file
        double FFTTime = i * (double)(1.0 / samplerateAudio);
        complex<double> power(0, (double)(-2.0 * PI * FFTTime * frequency));
        double w = min(2.0 * i / sampleCount, 2.0 - (2.0 * i / sampleCount));
        //calc += abs(PCM[time + i] * pow(E, power));
        calc += w * PCM[time + i] * pow(E, power);
    }
    
    return abs(calc) / (double)cnt;
}
void makeTxt(vector<vector<vector<unsigned char>>> & v, int x, int y, int multi, int type)
{
    for(int i = 0; i < 3; i++)
        for(int j = 0; j < 5; j++)
        {
            int ex = i * multi + x, exe = i * multi + x + multi;
            int ey = y + (4 - j) * multi, eye = y + (4 - j) * multi + multi;
            
            for(int k = ex; k < exe; k++)
                for(register int l = ey; l < eye; l++)
                    for(register int a = 0; a < 3; a++)
                    {
                        if(txtDot[type][j][i] == '#')
                            v[k][l][a] = 0;
                        else
                            v[k][l][a] = 255;
                    }
        }
}
void makeNum(vector<vector<vector<unsigned char>>> & v, int x, int y, int multi, int type)
{
    for(int i = 0; i < 3; i++)
        for(int j = 0; j < 5; j++)
        {
            int ex = i * multi + x, exe = i * multi + x + multi;
            int ey = y + (4 - j) * multi, eye = y + (4 - j) * multi + multi;
            
            for(int k = ex; k < exe; k++)
                for(register int l = ey; l < eye; l++)
                    for(register int a = 0; a < 3; a++)
                    {
                        if(numDot[type][j][i] == '#')
                            v[k][l][a] = 0;
                        else
                            v[k][l][a] = 255;
                    }
        }
}
void makeImage(vector<int> & dft)
{
    char direct[100];
    sprintf(direct, "outputimg/img%05d.bmp", dft.size());
    FILE *fp2 = fopen(direct, "wb");
    //vector<unsigned char> v2(3, 255);
    //vector<vector<unsigned char>> v1(720, v2);
    // B G R
    //vector<vector<vector<unsigned char>>> v(1280, v1);
    vector<vector<vector<unsigned char>>> v = imageOriginal;
    
    string code = pitchTable[dft.back() % 12];
    int pitch = dft.back() / 12;
    //cout << code << pitch << endl;
    
    //code and pitch making
    //(20, 20) ~ (79, 119), (100, 20) ~ (159, 119), (180, 20) ~ (239, 119), (860, 20) ~ (919, 119)
    if(code.size() == 2)
    {
        makeTxt(v, 120, 20, 20, code[0] - 'A');
        makeTxt(v, 200, 20, 20, 7);
    }
    else
        makeTxt(v, 120, 20, 20, code[0] - 'A');
    makeNum(v, 20, 20, 20, pitch);
    // 0: jinseong, 1: amaehada, 2: gaseong
    
    // graph teduri making
    //(20, 680) ~ (39, 699), (20, 180) ~ (39, 199), (1240, 680) ~ (1259, 699), (1240, 180) ~ (1259, 199)
    
    for(int i = 0; i < 60; i++)
    {
        int time = dft.size() - i - 1;
        if(0 > time)
            break;
        
        // x = 20px, y = 10px
        int ex = 1220 - 20 * i, ey = 200 + 10 * dft[time];
        for(int i = ex; i < ex + 20; i++)
            for(register int j = ey; j < ey + 10; j++)
                for(register int k = 0; k < 3; k++)
                    v[i][j][k] = v[1279][j][k];
    }
    
    for(register int i = 0; i < 54; i++)
        putc(imageBUFFER[i], fp2);
    
    for(int i = 0; i < 720; i++)
        for(register int j = 0; j < 1280; j++)
            for(register int k = 0; k < 3; k++)
                putc(v[j][i][k], fp2);

    fclose(fp2);
}

int main()
{
    auto startWav = chrono::high_resolution_clock::now();
    if(!system(NULL)) // requires FFmpeg
        exit(EXIT_FAILURE);
    
    // header reading
    FILE *fp2 = fopen("dummy/header.bmp", "rb");
    for(int i = 0; i < 54; i++)
        imageBUFFER[i] = getc(fp2);
    
    // original image vector making
    for(int i = 0; i < 1280; i++)
    {
        vector<unsigned char> v2(3, 0);
        vector<vector<unsigned char>> v1(720, v2);
        imageOriginal.push_back(v1);
    }
    
    // original image vector reading
    for(int i = 0; i < 720; i++)
        for(int j = 0; j < 1280; j++)
            for(int k = 0; k < 3; k++)
                imageOriginal[j][i][k] = getc(fp2);
    
    fclose(fp2);
    
    // input freq range
    cout << "Input min Freq and Max Freq (table: in freqtable.txt)" << endl;
    cin >> minFreq >> maxFreq;
    
    // header setting
    // change fileName in line 279, 354 / frameRate 318, 339, 353
    system("ffmpeg -y -i 440Hz.mp3 dummy\\audio.wav");
    FILE *fp1 = fopen("dummy/audio.wav", "rb");
    WAV head;
    WAVt headt;
    fread(&head, sizeof(head), 1, fp1);
    headt = readheader(head);
    samplerateAudio = headt.SampleRate;
    channelsAudio = headt.NumChannels;
    bytesAudio = headt.BlockAlign / channelsAudio;
    unsigned char BUFFER[4];
    vector<double> PCM;
    vector<int> dft, freq(48, 48);
    //cout << headt.NumChannels << " " << headt.BlockAlign << " " << headt.BitsPerSample  << endl;
    
    // Audio input
    while(1)
    {
        if(feof(fp1))
            break;
        
        double s = 0;
        for(int j = 0; j < channelsAudio; j++)
        {
            fread(BUFFER, bytesAudio, 1, fp1);
            s += ler(BUFFER, bytesAudio);
        }
        PCM.push_back(s);
    }
    lengthAudio = PCM.size();
    fclose(fp1);
    
    cout << "File I/O end" << endl;
    
    // calc freq (
    vector<double> freqency(48);
    for(int i = 0; i < 48; i++)
        freqency[i] = 440.0 * pow(pow(2, 1 /(double)12), i - 33);

    // do dft in every 1/frame s, freq
    for(int i = 0; i < lengthAudio; i += (samplerateAudio / 10))
    {
        vector<double> dftResult(48);
        pair<int, double> dftMax = {0, 0};
        for(int j = minFreq; j <= maxFreq; j++)
        {
            dftResult[j] = calcFourier(PCM, i, freqency[j]);
            if(dftMax.second <= dftResult[j])
                dftMax = {j, dftResult[j]};
        }
        if(33 <= dftMax.first)
        {
            double freqd2 = dftResult[dftMax.first - 12];
            if(freqd2 >= dftMax.second * 0.8)
                dftMax = {dftMax.first - 12, freqd2};
        }

        dft.push_back(dftMax.first);
        makeImage(dft);
        
        // see percentage
        if(i % samplerateAudio < samplerateAudio / 10)
        {
            auto endWav = chrono::high_resolution_clock::now();
            auto duration = chrono::duration_cast<chrono::milliseconds>(endWav - startWav).count();
            cout << "now making " << i / samplerateAudio << "s" << endl << "Rendering took " << duration << " ms." <<
            endl;
        }
    }
    
    auto endWav = chrono::high_resolution_clock::now();
    auto duration = chrono::duration_cast<chrono::milliseconds>(endWav - startWav).count();
    cout << endl << "Rendering done " << duration << " ms." << endl << "Start making video" << endl;
    
    // re encoding
    system("ffmpeg -y -r 10 -i outputimg\\img%05d.bmp dummy\\output.mkv");
    system("ffmpeg -y -i dummy\\output.mkv -i 440Hz.mp3 -c:v copy -c:a copy dummy\\finaloutput.mp4");
    system("ffmpeg -y -i dummy\\finaloutput.mp4 -vcodec libx264 -b 3000k -c:a copy finaloutput.mp4");
    
    endWav = chrono::high_resolution_clock::now();
    duration = chrono::duration_cast<chrono::milliseconds>(endWav - startWav).count();
    cout << endl << "Rendering done " << duration << " ms." << endl;
    cout << "Please erase images in outputimg folder" << endl;
    
    return 0;
}