//
// Created by andgra on 13.10.2018.
//

#ifndef STRAIGHT_RAY_EXO_H
#define STRAIGHT_RAY_EXO_H

#include "comp.h"

// Decompiled with JetBrains decompiler
// Type: Exocortex.DSP.Fourier
// Assembly: Exocortex.DSP.v1, Version=1.0.1373.39154, Culture=neutral, PublicKeyToken=null
// MVID: 5CC5C8AF-4461-49AD-8B77-23F6A0FE5383
// Assembly location: D:\универ\Магистерская\TracingUltrasound\TracingUltrasound 07.10.18\TracingUltrasound\Exocortex.DSP.v1.dll


class Fourier {
private:
    static int** _reversedBits = new int*[12];
    static int** _reverseBits;
    static int _lookupTabletLength = -1;
    static double[,][]
    _uRLookup = (double[,
    ][])
    null;
    static double[,][]
    _uILookup = (double[,
    ][])
    null;
    static float[,][]
    _uRLookupF = (float[,
    ][])
    null;
    static float[,][]
    _uILookupF = (float[,
    ][])
    null;
    static bool _bufferFLocked = false;
    static float[]
    _bufferF = new float[0];
    static bool _bufferCFLocked = false;
    static compF[]
    _bufferCF = new compF[0];
    static bool _bufferCLocked = false;
    static comp *_bufferC = new comp[0];
    const int cMaxLength = 4096;
    const int cMinLength = 1;
    const int cMaxBits = 12;
    const int cMinBits = 0;

    Fourier() {
    }

    static void Swap(float& a, float& b) {
        float num = a;
        a = b;
        b = num;
    }

    static void Swap(double& a, double b) {
        double num = a;
        a = b;
        b = num;
    }

    static void Swap(compF& a, compF& b) {
        compF compF = a;
        a = b;
        b = compF;
    }

    static void Swap(comp& a, comp& b) {
        comp comp = a;
        a = b;
        b = comp;
    }

    static bool IsPowerOf2(int x) {
        return (x & (x - 1)) == 0;
    }

    static int Pow2(int exponent) {
        if (exponent >= 0 && exponent < 31)
            return 1 << exponent;
        return 0;
    }

    static int Log2(int x) {
        if (x <= 65536) {
            if (x <= 256) {
                if (x <= 16) {
                    if (x <= 4) {
                        if (x > 2)
                            return 2;
                        return x <= 1 ? 0 : 1;
                    }
                    return x <= 8 ? 3 : 4;
                }
                if (x <= 64)
                    return x <= 32 ? 5 : 6;
                return x <= 128 ? 7 : 8;
            }
            if (x <= 4096) {
                if (x <= 1024)
                    return x <= 512 ? 9 : 10;
                return x <= 2048 ? 11 : 12;
            }
            if (x <= 16384)
                return x <= 8192 ? 13 : 14;
            return x <= 32768 ? 15 : 16;
        }
        if (x <= 16777216) {
            if (x <= 1048576) {
                if (x <= 262144)
                    return x <= 131072 ? 17 : 18;
                return x <= 524288 ? 19 : 20;
            }
            if (x <= 4194304)
                return x <= 2097152 ? 21 : 22;
            return x <= 8388608 ? 23 : 24;
        }
        if (x <= 268435456) {
            if (x <= 67108864)
                return x <= 33554432 ? 25 : 26;
            return x <= 134217728 ? 27 : 28;
        }
        if (x > 1073741824)
            return 31;
        return x <= 536870912 ? 29 : 30;
    }

    static int ReverseBits(int index, int numberOfBits) {
        int num1 = 0;
        int num2 = 0;
        while (num2 < numberOfBits) {
            num1 = num1 << 1 | index & 1;
            index >>= 1;
            ++num2;
        }
        return num1;
    }

    static int[]

    GetReversedBits(int numberOfBits) {
        if (_reversedBits[numberOfBits - 1] == null) {
            int length = Pow2(numberOfBits);
            int[]
            numArray = new int[length];
            int index = 0;
            while (index < length) {
                int num1 = index;
                int num2 = 0;
                int num3 = 0;
                while (num3 < numberOfBits) {
                    num2 = num2 << 1 | num1 & 1;
                    num1 >>= 1;
                    checked{++num3;}
                }
                numArray[index] = num2;
                checked{++index;}
            }
            _reversedBits[checked(numberOfBits - 1)] = numArray;
        }
        return _reversedBits[checked(numberOfBits - 1)];
    }

    static void ReorderArray(float[] data) {
        Debug.Assert(data != null);
        int x = data.Length / 2;
        Debug.Assert(IsPowerOf2(x));
        Debug.Assert(x >= 1);
        Debug.Assert(x <= 4096);
        int[]
        reversedBits = GetReversedBits(Log2(x));
        int index = 0;
        while (index < x) {
            int num = reversedBits[index];
            if (num > index) {
                Swap(ref
                data[index << 1], ref
                data[num << 1]);
                Swap(ref
                data[checked(index << 1 + 1)], ref
                data[checked(num << 1 + 1)]);
            }
            checked{++index;}
        }
    }

    static void ReorderArray(double[] data) {
        Debug.Assert(data != null);
        int x = data.Length / 2;
        Debug.Assert(IsPowerOf2(x));
        Debug.Assert(x >= 1);
        Debug.Assert(x <= 4096);
        int[]
        reversedBits = GetReversedBits(Log2(x));
        int index = 0;
        while (index < x) {
            int num = reversedBits[index];
            if (num > index) {
                Swap(ref
                data[index << 1], ref
                data[num << 1]);
                Swap(ref
                data[index << 2], ref
                data[num << 2]);
            }
            checked{++index;}
        }
    }

    static void ReorderArray(comp *data) {
        Debug.Assert(data != null);
        int length = data.Length;
        Debug.Assert(IsPowerOf2(length));
        Debug.Assert(length >= 1);
        Debug.Assert(length <= 4096);
        int[]
        reversedBits = GetReversedBits(Log2(length));
        int index1 = 0;
        while (index1 < length) {
            int index2 = reversedBits[index1];
            if (index2 > index1) {
                comp comp = data[index1];
                data[index1] = data[index2];
                data[index2] = comp;
            }
            checked{++index1;}
        }
    }

    static void ReorderArray(compF[] data) {
        Debug.Assert(data != null);
        int length = data.Length;
        Debug.Assert(IsPowerOf2(length));
        Debug.Assert(length >= 1);
        Debug.Assert(length <= 4096);
        int[]
        reversedBits = GetReversedBits(Log2(length));
        int index1 = 0;
        while (index1 < length) {
            int index2 = reversedBits[index1];
            if (index2 > index1) {
                compF compF = data[index1];
                data[index1] = data[index2];
                data[index2] = compF;
            }
            checked{++index1;}
        }
    }

    static int _ReverseBits(int bits, int n) {
        int num1 = 0;
        int num2 = 0;
        while (num2 < n) {
            num1 = num1 << 1 | bits & 1;
            bits >>= 1;
            checked{++num2;}
        }
        return num1;
    }

    static void InitializeReverseBits(int levels) {
        _reverseBits = new int[checked(levels + 1)][];
        int n = 0;
        while (n < checked(levels + 1)) {
            int length = checked((int) Math.Pow(2.0, unchecked((double) n)));
            _reverseBits[n] = new int[length];
            int bits = 0;
            while (bits < length) {
                _reverseBits[n][bits] = _ReverseBits(bits, n);
                checked{++bits;}
            }
            checked{++n;}
        }
    }

    static void SyncLookupTableLength(int length) {
        Debug.Assert(length < 10240);
        Debug.Assert(length >= 0);
        if (length <= _lookupTabletLength)
            return;
        int levels = checked((int) Math.Ceiling(Math.Log(unchecked((double) length), 2.0)));
        InitializeReverseBits(levels);
        InitializecompRotations(levels);
        _lookupTabletLength = length;
    }

    static int GetLookupTableLength() {
        return _lookupTabletLength;
    }

    static void ClearLookupTables() {
        _uRLookup = (double[,][]) null;
        _uILookup = (double[,][]) null;
        _uRLookupF = (float[,][]) null;
        _uILookupF = (float[,][]) null;
        _lookupTabletLength = -1;
    }

    static void InitializecompRotations(int levels) {
        int num1 = levels;
        _uRLookup = new double[checked(levels + 1), 2][];
        _uILookup = new double[checked(levels + 1), 2][];
        _uRLookupF = new float[checked(levels + 1), 2][];
        _uILookupF = new float[checked(levels + 1), 2][];
        int num2 = 1;
        int index1 = 1;
        while (index1 <= num1) {
            int length = num2;
            num2 <<= 1;
            double num3 = 1.0;
            double num4 = 0.0;
            double num5 = Math.PI / (double) length * 1.0;
            double num6 = Math.Cos(num5);
            double num7 = Math.Sin(num5);
            _uRLookup[index1, 0] = new double[length];
            _uILookup[index1, 0] = new double[length];
            _uRLookupF[index1, 0] = new float[length];
            _uILookupF[index1, 0] = new float[length];
            int index2 = 0;
            while (index2 < length) {
                _uRLookupF[index1, 0][index2] = (float) (_uRLookup[index1, 0][index2] = num3);
                _uILookupF[index1, 0][index2] = (float) (_uILookup[index1, 0][index2] = num4);
                double num8 = num3 * num7 + num4 * num6;
                num3 = num3 * num6 - num4 * num7;
                num4 = num8;
                checked{++index2;}
            }
            double num9 = 1.0;
            double num10 = 0.0;
            double num11 = Math.PI / (double) length * -1.0;
            double num12 = Math.Cos(num11);
            double num13 = Math.Sin(num11);
            _uRLookup[index1, 1] = new double[length];
            _uILookup[index1, 1] = new double[length];
            _uRLookupF[index1, 1] = new float[length];
            _uILookupF[index1, 1] = new float[length];
            int index3 = 0;
            while (index3 < length) {
                _uRLookupF[index1, 1][index3] = (float) (_uRLookup[index1, 1][index3] = num9);
                _uILookupF[index1, 1][index3] = (float) (_uILookup[index1, 1][index3] = num10);
                double num8 = num9 * num13 + num10 * num12;
                num9 = num9 * num12 - num10 * num13;
                num10 = num8;
                checked{++index3;}
            }
            checked{++index1;}
        }
    }

    static void LockBufferF(int length, ref float[] buffer) {
        Debug.Assert(!_bufferFLocked);
        _bufferFLocked = true;
        if (length >= _bufferF.Length)
            _bufferF = new float[length];
        buffer = _bufferF;
    }

    static void UnlockBufferF(ref float[] buffer) {
        Debug.Assert(_bufferF == buffer);
        Debug.Assert(_bufferFLocked);
        _bufferFLocked = false;
        buffer = (float[]) null;
    }

    static void LinearFFT(float[] data, int start, int inc, int length, FourierDirection direction) {
        Debug.Assert(data != null);
        Debug.Assert(start >= 0);
        Debug.Assert(inc >= 1);
        Debug.Assert(length >= 1);
        Debug.Assert(checked(start + inc * (length - 1) * 2) < data.Length);
        float[]
        buffer = (float[]) null;
        LockBufferF(checked(length * 2), ref
        buffer);
        int index1 = start;
        int index2 = 0;
        while (index2 < checked(length * 2)) {
            buffer[index2] = data[index1];
            checked{index1 += inc;}
            checked{++index2;}
        }
        FFT(buffer, length, direction);
        int index3 = start;
        int index4 = 0;
        while (index4 < length) {
            data[index3] = buffer[index4];
            checked{index3 += inc;}
            checked{++index4;}
        }
        UnlockBufferF(ref
        buffer);
    }

    static void LinearFFT_Quick(float[] data, int start, int inc, int length, FourierDirection direction) {
        float[]
        buffer = (float[]) null;
        LockBufferF(checked(length * 2), ref
        buffer);
        int index1 = start;
        int index2 = 0;
        while (index2 < checked(length * 2)) {
            buffer[index2] = data[index1];
            checked{index1 += inc;}
            checked{++index2;}
        }
        FFT_Quick(buffer, length, direction);
        int index3 = start;
        int index4 = 0;
        while (index4 < length) {
            data[index3] = buffer[index4];
            checked{index3 += inc;}
            checked{++index4;}
        }
        UnlockBufferF(ref
        buffer);
    }

    static void LockBufferCF(int length, ref compF[] buffer) {
        Debug.Assert(length >= 0);
        Debug.Assert(!_bufferCFLocked);
        _bufferCFLocked = true;
        if (length != _bufferCF.Length)
            _bufferCF = new compF[length];
        buffer = _bufferCF;
    }

    static void UnlockBufferCF(ref compF[] buffer) {
        Debug.Assert(_bufferCF == buffer);
        Debug.Assert(_bufferCFLocked);
        _bufferCFLocked = false;
        buffer = (compF[])
        null;
    }

    static void LinearFFT(compF[] data, int start, int inc, int length, FourierDirection direction) {
        Debug.Assert(data != null);
        Debug.Assert(start >= 0);
        Debug.Assert(inc >= 1);
        Debug.Assert(length >= 1);
        Debug.Assert(checked(start + inc * (length - 1)) < data.Length);
        compF[]
        buffer = (compF[]) null;
        LockBufferCF(length, ref
        buffer);
        int index1 = start;
        int index2 = 0;
        while (index2 < length) {
            buffer[index2] = data[index1];
            checked{index1 += inc;}
            checked{++index2;}
        }
        FFT(buffer, length, direction);
        int index3 = start;
        int index4 = 0;
        while (index4 < length) {
            data[index3] = buffer[index4];
            checked{index3 += inc;}
            checked{++index4;}
        }
        UnlockBufferCF(ref
        buffer);
    }

    static void LinearFFT_Quick(compF[] data, int start, int inc, int length, FourierDirection direction) {
        compF[]
        buffer = (compF[]) null;
        LockBufferCF(length, ref
        buffer);
        int index1 = start;
        int index2 = 0;
        while (index2 < length) {
            buffer[index2] = data[index1];
            checked{index1 += inc;}
            checked{++index2;}
        }
        FFT(buffer, length, direction);
        int index3 = start;
        int index4 = 0;
        while (index4 < length) {
            data[index3] = buffer[index4];
            checked{index3 += inc;}
            checked{++index4;}
        }
        UnlockBufferCF(ref
        buffer);
    }

    static void LockBufferC(int length, ref comp

    * buffer)
    {
        Debug.Assert(length >= 0);
        Debug.Assert(!_bufferCLocked);
        _bufferCLocked = true;
        if (length >= _bufferC.Length)
            _bufferC = new comp[length];
        buffer = _bufferC;
    }

    static void UnlockBufferC(ref comp

    * buffer)
    {
        Debug.Assert(_bufferC == buffer);
        Debug.Assert(_bufferCLocked);
        _bufferCLocked = false;
        buffer = (comp *) null;
    }

    static void LinearFFT(comp *data, int start, int inc, int length, FourierDirection direction) {
        Debug.Assert(data != null);
        Debug.Assert(start >= 0);
        Debug.Assert(inc >= 1);
        Debug.Assert(length >= 1);
        Debug.Assert(checked(start + inc * (length - 1)) < data.Length);
        comp *buffer = (comp *) null;
        LockBufferC(length, ref
        buffer);
        int index1 = start;
        int index2 = 0;
        while (index2 < length) {
            buffer[index2] = data[index1];
            checked{index1 += inc;}
            checked{++index2;}
        }
        FFT(buffer, length, direction);
        int index3 = start;
        int index4 = 0;
        while (index4 < length) {
            data[index3] = buffer[index4];
            checked{index3 += inc;}
            checked{++index4;}
        }
        UnlockBufferC(ref
        buffer);
    }

    static void LinearFFT_Quick(comp *data, int start, int inc, int length, FourierDirection direction) {
        comp *buffer = (comp *) null;
        LockBufferC(length, ref
        buffer);
        int index1 = start;
        int index2 = 0;
        while (index2 < length) {
            buffer[index2] = data[index1];
            checked{index1 += inc;}
            checked{++index2;}
        }
        FFT_Quick(buffer, length, direction);
        int index3 = start;
        int index4 = 0;
        while (index4 < length) {
            data[index3] = buffer[index4];
            checked{index3 += inc;}
            checked{++index4;}
        }
        UnlockBufferC(ref
        buffer);
    }

public:
    static void FFT(float[] data, int length, FourierDirection direction) {
        Debug.Assert(data != null);
        Debug.Assert(data.Length >= checked(length * 2));
        Debug.Assert(IsPowerOf2(length));
        SyncLookupTableLength(length);
        int num1 = Log2(length);
        ReorderArray(data);
        int num2 = 1;
        int index1 = direction == FourierDirection.Forward ? 0 : 1;
        int index2 = 1;
        while (index2 <= num1) {
            int num3 = num2;
            num2 <<= 1;
            float[]
            numArray1 = _uRLookupF[index2, index1];
            float[]
            numArray2 = _uILookupF[index2, index1];
            int index3 = 0;
            while (index3 < num3) {
                float num4 = numArray1[index3];
                float num5 = numArray2[index3];
                int num6 = index3;
                while (num6 < length) {
                    int index4 = num6 << 1;
                    int index5 = checked(num6 + num3) << 1;
                    float num7 = data[index5];
                    float num8 = data[checked(index5 + 1)];
                    float num9 = (float) ((double) num7 * (double) num4 - (double) num8 * (double) num5);
                    float num10 = (float) ((double) num7 * (double) num5 + (double) num8 * (double) num4);
                    float num11 = data[index4];
                    float num12 = data[checked(index4 + 1)];
                    data[index4] = num11 + num9;
                    data[checked(index4 + 1)] = num12 + num10;
                    data[index5] = num11 - num9;
                    data[checked(index5 + 1)] = num12 - num10;
                    checked{num6 += num2;}
                }
                checked{++index3;}
            }
            checked{++index2;}
        }
    }

    static void FFT_Quick(float[] data, int length, FourierDirection direction) {
        int num1 = Log2(length);
        ReorderArray(data);
        int num2 = 1;
        int index1 = direction == FourierDirection.Forward ? 0 : 1;
        int index2 = 1;
        while (index2 <= num1) {
            int num3 = num2;
            num2 <<= 1;
            float[]
            numArray1 = _uRLookupF[index2, index1];
            float[]
            numArray2 = _uILookupF[index2, index1];
            int index3 = 0;
            while (index3 < num3) {
                float num4 = numArray1[index3];
                float num5 = numArray2[index3];
                int num6 = index3;
                while (num6 < length) {
                    int index4 = num6 << 1;
                    int index5 = checked(num6 + num3) << 1;
                    float num7 = data[index5];
                    float num8 = data[checked(index5 + 1)];
                    float num9 = (float) ((double) num7 * (double) num4 - (double) num8 * (double) num5);
                    float num10 = (float) ((double) num7 * (double) num5 + (double) num8 * (double) num4);
                    float num11 = data[index4];
                    float num12 = data[checked(index4 + 1)];
                    data[index4] = num11 + num9;
                    data[checked(index4 + 1)] = num12 + num10;
                    data[index5] = num11 - num9;
                    data[checked(index5 + 1)] = num12 - num10;
                    checked{num6 += num2;}
                }
                checked{++index3;}
            }
            checked{++index2;}
        }
    }

    static void FFT(compF[] data, int length, FourierDirection direction) {
        if (data == null)
            throw new ArgumentNullException(nameof(data));
        if (data.Length < length)
            throw new ArgumentOutOfRangeException(nameof(length), (object) length, "must be at least as large as 'data.Length' parameter");
        if (!IsPowerOf2(length))
            throw new ArgumentOutOfRangeException(nameof(length), (object) length, "must be a power of 2");
        SyncLookupTableLength(length);
        int num1 = Log2(length);
        ReorderArray(data);
        int num2 = 1;
        int index1 = direction == FourierDirection.Forward ? 0 : 1;
        int index2 = 1;
        while (index2 <= num1) {
            int num3 = num2;
            num2 <<= 1;
            float[]
            numArray1 = _uRLookupF[index2, index1];
            float[]
            numArray2 = _uILookupF[index2, index1];
            int index3 = 0;
            while (index3 < num3) {
                float num4 = numArray1[index3];
                float num5 = numArray2[index3];
                int index4 = index3;
                while (index4 < length) {
                    int index5 = checked(index4 + num3);
                    float re1 = data[index5].Re;
                    float im1 = data[index5].Im;
                    float num6 = (float) ((double) re1 * (double) num4 - (double) im1 * (double) num5);
                    float num7 = (float) ((double) re1 * (double) num5 + (double) im1 * (double) num4);
                    float re2 = data[index4].Re;
                    float im2 = data[index4].Im;
                    data[index4].Re = re2 + num6;
                    data[index4].Im = im2 + num7;
                    data[index5].Re = re2 - num6;
                    data[index5].Im = im2 - num7;
                    checked{index4 += num2;}
                }
                checked{++index3;}
            }
            checked{++index2;}
        }
    }

    static void FFT_Quick(compF[] data, int length, FourierDirection direction) {
        int num1 = Log2(length);
        ReorderArray(data);
        int num2 = 1;
        int index1 = direction == FourierDirection.Forward ? 0 : 1;
        int index2 = 1;
        while (index2 <= num1) {
            int num3 = num2;
            num2 <<= 1;
            float[]
            numArray1 = _uRLookupF[index2, index1];
            float[]
            numArray2 = _uILookupF[index2, index1];
            int index3 = 0;
            while (index3 < num3) {
                float num4 = numArray1[index3];
                float num5 = numArray2[index3];
                int index4 = index3;
                while (index4 < length) {
                    int index5 = checked(index4 + num3);
                    float re1 = data[index5].Re;
                    float im1 = data[index5].Im;
                    float num6 = (float) ((double) re1 * (double) num4 - (double) im1 * (double) num5);
                    float num7 = (float) ((double) re1 * (double) num5 + (double) im1 * (double) num4);
                    float re2 = data[index4].Re;
                    float im2 = data[index4].Im;
                    data[index4].Re = re2 + num6;
                    data[index4].Im = im2 + num7;
                    data[index5].Re = re2 - num6;
                    data[index5].Im = im2 - num7;
                    checked{index4 += num2;}
                }
                checked{++index3;}
            }
            checked{++index2;}
        }
    }

    static void FFT(compF[] data, FourierDirection direction) {
        if (data == null)
            throw new ArgumentNullException(nameof(data));
        FFT(data, data.Length, direction);
    }

    static void FFT(comp *data, int length, FourierDirection direction) {
        if (data == null)
            throw new ArgumentNullException(nameof(data));
        if (data.Length < length)
            throw new ArgumentOutOfRangeException(nameof(length), (object) length, "must be at least as large as 'data.Length' parameter");
        if (!IsPowerOf2(length))
            throw new ArgumentOutOfRangeException(nameof(length), (object) length, "must be a power of 2");
        SyncLookupTableLength(length);
        int num1 = Log2(length);
        ReorderArray(data);
        int num2 = 1;
        int index1 = direction == FourierDirection.Forward ? 0 : 1;
        int index2 = 1;
        while (index2 <= num1) {
            int num3 = num2;
            num2 <<= 1;
            double[]
            numArray1 = _uRLookup[index2, index1];
            double[]
            numArray2 = _uILookup[index2, index1];
            int index3 = 0;
            while (index3 < num3) {
                double num4 = numArray1[index3];
                double num5 = numArray2[index3];
                int index4 = index3;
                while (index4 < length) {
                    int index5 = checked(index4 + num3);
                    double re1 = data[index5].Re;
                    double im1 = data[index5].Im;
                    double num6 = re1 * num4 - im1 * num5;
                    double num7 = re1 * num5 + im1 * num4;
                    double re2 = data[index4].Re;
                    double im2 = data[index4].Im;
                    data[index4].Re = re2 + num6;
                    data[index4].Im = im2 + num7;
                    data[index5].Re = re2 - num6;
                    data[index5].Im = im2 - num7;
                    checked{index4 += num2;}
                }
                checked{++index3;}
            }
            checked{++index2;}
        }
    }

    static void FFT_Quick(comp *data, int length, FourierDirection direction) {
        int num1 = Log2(length);
        ReorderArray(data);
        int num2 = 1;
        int index1 = direction == FourierDirection.Forward ? 0 : 1;
        int index2 = 1;
        while (index2 <= num1) {
            int num3 = num2;
            num2 <<= 1;
            double[]
            numArray1 = _uRLookup[index2, index1];
            double[]
            numArray2 = _uILookup[index2, index1];
            int index3 = 0;
            while (index3 < num3) {
                double num4 = numArray1[index3];
                double num5 = numArray2[index3];
                int index4 = index3;
                while (index4 < length) {
                    int index5 = checked(index4 + num3);
                    double re1 = data[index5].Re;
                    double im1 = data[index5].Im;
                    double num6 = re1 * num4 - im1 * num5;
                    double num7 = re1 * num5 + im1 * num4;
                    double re2 = data[index4].Re;
                    double im2 = data[index4].Im;
                    data[index4].Re = re2 + num6;
                    data[index4].Im = im2 + num7;
                    data[index5].Re = re2 - num6;
                    data[index5].Im = im2 - num7;
                    checked{index4 += num2;}
                }
                checked{++index3;}
            }
            checked{++index2;}
        }
    }

    static void RFFT(float[] data, FourierDirection direction) {
        if (data == null)
            throw new ArgumentNullException(nameof(data));
        RFFT(data, data.Length, direction);
    }

    static void RFFT(float[] data, int length, FourierDirection direction) {
        if (data == null)
            throw new ArgumentNullException(nameof(data));
        if (data.Length < length)
            throw new ArgumentOutOfRangeException(nameof(length), (object) length, "must be at least as large as 'data.Length' parameter");
        if (!IsPowerOf2(length))
            throw new ArgumentOutOfRangeException(nameof(length), (object) length, "must be a power of 2");
        float num1 = 0.5f;
        float num2 = 3.141593f / (float) (length / 2);
        float num3;
        if (direction == FourierDirection.Forward) {
            num3 = -0.5f;
            FFT(data, length / 2, direction);
        } else {
            num3 = 0.5f;
            num2 = -num2;
        }
        float num4 = (float) Math.Sin(0.5 * (double) num2);
        float num5 = -2f * num4 * num4;
        float num6 = (float) Math.Sin((double) num2);
        float num7 = 1f + num5;
        float num8 = num6;
        int num9 = 1;
        while (num9 < length / 4) {
            int index1 = checked(2 * num9);
            int index2 = checked(length - 2 * num9);
            float num10 = num1 * (data[index1] + data[index2]);
            float num11 = num1 * (data[checked(index1 + 1)] - data[checked(index2 + 1)]);
            float num12 = (float) (-(double) num3 * ((double) data[checked(index1 + 1)] + (double) data[checked(index2 + 1)]));
            float num13 = num3 * (data[index1] - data[index2]);
            data[index1] = (float) ((double) num10 + (double) num7 * (double) num12 - (double) num8 * (double) num13);
            data[checked(index1 + 1)] = (float) ((double) num11 + (double) num7 * (double) num13 + (double) num8 * (double) num12);
            data[index2] = (float) ((double) num10 - (double) num7 * (double) num12 + (double) num8 * (double) num13);
            data[checked(index2 + 1)] = (float) (-(double) num11 + (double) num7 * (double) num13 + (double) num8 * (double) num12);
            float num14;
            num7 = (float) ((double) (num14 = num7) * (double) num5 - (double) num8 * (double) num6) + num7;
            num8 = (float) ((double) num8 * (double) num5 + (double) num14 * (double) num6) + num8;
            checked{++num9;}
        }
        if (direction == FourierDirection.Forward) {
            float num10 = data[0];
            data[0] = num10 + data[1];
            data[1] = num10 - data[1];
        } else {
            float num10 = data[0];
            data[0] = num1 * (num10 + data[1]);
            data[1] = num1 * (num10 - data[1]);
            FFT(data, length / 2, direction);
        }
    }

    static void FFT2(float[] data, int xLength, int yLength, FourierDirection direction) {
        if (data == null)
            throw new ArgumentNullException(nameof(data));
        if (data.Length < checked(xLength * yLength * 2))
            throw new ArgumentOutOfRangeException("data.Length", (object) data.Length, "must be at least as large as 'xLength * yLength * 2' parameter");
        if (!IsPowerOf2(xLength))
            throw new ArgumentOutOfRangeException(nameof(xLength), (object) xLength, "must be a power of 2");
        if (!IsPowerOf2(yLength))
            throw new ArgumentOutOfRangeException(nameof(yLength), (object) yLength, "must be a power of 2");
        int inc1 = 1;
        int inc2 = xLength;
        if (xLength > 1) {
            SyncLookupTableLength(xLength);
            int num = 0;
            while (num < yLength) {
                int start = checked(num * inc2);
                LinearFFT_Quick(data, start, inc1, xLength, direction);
                checked{++num;}
            }
        }
        if (yLength <= 1)
            return;
        SyncLookupTableLength(yLength);
        int num1 = 0;
        while (num1 < xLength) {
            int start = checked(num1 * inc1);
            LinearFFT_Quick(data, start, inc2, yLength, direction);
            checked{++num1;}
        }
    }

    static void FFT2(compF[] data, int xLength, int yLength, FourierDirection direction) {
        if (data == null)
            throw new ArgumentNullException(nameof(data));
        if (data.Length < checked(xLength * yLength))
            throw new ArgumentOutOfRangeException("data.Length", (object) data.Length, "must be at least as large as 'xLength * yLength' parameter");
        if (!IsPowerOf2(xLength))
            throw new ArgumentOutOfRangeException(nameof(xLength), (object) xLength, "must be a power of 2");
        if (!IsPowerOf2(yLength))
            throw new ArgumentOutOfRangeException(nameof(yLength), (object) yLength, "must be a power of 2");
        int inc1 = 1;
        int inc2 = xLength;
        if (xLength > 1) {
            SyncLookupTableLength(xLength);
            int num = 0;
            while (num < yLength) {
                int start = checked(num * inc2);
                LinearFFT_Quick(data, start, inc1, xLength, direction);
                checked{++num;}
            }
        }
        if (yLength <= 1)
            return;
        SyncLookupTableLength(yLength);
        int num1 = 0;
        while (num1 < xLength) {
            int start = checked(num1 * inc1);
            LinearFFT_Quick(data, start, inc2, yLength, direction);
            checked{++num1;}
        }
    }

    static void FFT2(comp *data, int xLength, int yLength, FourierDirection direction) {
        if (data == null)
            throw new ArgumentNullException(nameof(data));
        if (data.Length < checked(xLength * yLength))
            throw new ArgumentOutOfRangeException("data.Length", (object) data.Length, "must be at least as large as 'xLength * yLength' parameter");
        if (!IsPowerOf2(xLength))
            throw new ArgumentOutOfRangeException(nameof(xLength), (object) xLength, "must be a power of 2");
        if (!IsPowerOf2(yLength))
            throw new ArgumentOutOfRangeException(nameof(yLength), (object) yLength, "must be a power of 2");
        int inc1 = 1;
        int inc2 = xLength;
        if (xLength > 1) {
            SyncLookupTableLength(xLength);
            int num = 0;
            while (num < yLength) {
                int start = checked(num * inc2);
                LinearFFT_Quick(data, start, inc1, xLength, direction);
                checked{++num;}
            }
        }
        if (yLength <= 1)
            return;
        SyncLookupTableLength(yLength);
        int num1 = 0;
        while (num1 < xLength) {
            int start = checked(num1 * inc1);
            LinearFFT_Quick(data, start, inc2, yLength, direction);
            checked{++num1;}
        }
    }

    static void FFT3(compF[] data, int xLength, int yLength, int zLength, FourierDirection direction) {
        if (data == null)
            throw new ArgumentNullException(nameof(data));
        if (data.Length < checked(xLength * yLength * zLength))
            throw new ArgumentOutOfRangeException("data.Length", (object) data.Length, "must be at least as large as 'xLength * yLength * zLength' parameter");
        if (!IsPowerOf2(xLength))
            throw new ArgumentOutOfRangeException(nameof(xLength), (object) xLength, "must be a power of 2");
        if (!IsPowerOf2(yLength))
            throw new ArgumentOutOfRangeException(nameof(yLength), (object) yLength, "must be a power of 2");
        if (!IsPowerOf2(zLength))
            throw new ArgumentOutOfRangeException(nameof(zLength), (object) zLength, "must be a power of 2");
        int inc1 = 1;
        int inc2 = xLength;
        int inc3 = checked(xLength * yLength);
        if (xLength > 1) {
            SyncLookupTableLength(xLength);
            int num1 = 0;
            while (num1 < zLength) {
                int num2 = 0;
                while (num2 < yLength) {
                    int start = checked(num2 * inc2 + num1 * inc3);
                    LinearFFT_Quick(data, start, inc1, xLength, direction);
                    checked{++num2;}
                }
                checked{++num1;}
            }
        }
        if (yLength > 1) {
            SyncLookupTableLength(yLength);
            int num1 = 0;
            while (num1 < zLength) {
                int num2 = 0;
                while (num2 < xLength) {
                    int start = checked(num1 * inc3 + num2 * inc1);
                    LinearFFT_Quick(data, start, inc2, yLength, direction);
                    checked{++num2;}
                }
                checked{++num1;}
            }
        }
        if (zLength <= 1)
            return;
        SyncLookupTableLength(zLength);
        int num3 = 0;
        while (num3 < yLength) {
            int num1 = 0;
            while (num1 < xLength) {
                int start = checked(num3 * inc2 + num1 * inc1);
                LinearFFT_Quick(data, start, inc3, zLength, direction);
                checked{++num1;}
            }
            checked{++num3;}
        }
    }

    static void FFT3(comp *data, int xLength, int yLength, int zLength, FourierDirection direction) {
        if (data == null)
            throw new ArgumentNullException(nameof(data));
        if (data.Length < checked(xLength * yLength * zLength))
            throw new ArgumentOutOfRangeException("data.Length", (object) data.Length, "must be at least as large as 'xLength * yLength * zLength' parameter");
        if (!IsPowerOf2(xLength))
            throw new ArgumentOutOfRangeException(nameof(xLength), (object) xLength, "must be a power of 2");
        if (!IsPowerOf2(yLength))
            throw new ArgumentOutOfRangeException(nameof(yLength), (object) yLength, "must be a power of 2");
        if (!IsPowerOf2(zLength))
            throw new ArgumentOutOfRangeException(nameof(zLength), (object) zLength, "must be a power of 2");
        int inc1 = 1;
        int inc2 = xLength;
        int inc3 = checked(xLength * yLength);
        if (xLength > 1) {
            SyncLookupTableLength(xLength);
            int num1 = 0;
            while (num1 < zLength) {
                int num2 = 0;
                while (num2 < yLength) {
                    int start = checked(num2 * inc2 + num1 * inc3);
                    LinearFFT_Quick(data, start, inc1, xLength, direction);
                    checked{++num2;}
                }
                checked{++num1;}
            }
        }
        if (yLength > 1) {
            SyncLookupTableLength(yLength);
            int num1 = 0;
            while (num1 < zLength) {
                int num2 = 0;
                while (num2 < xLength) {
                    int start = checked(num1 * inc3 + num2 * inc1);
                    LinearFFT_Quick(data, start, inc2, yLength, direction);
                    checked{++num2;}
                }
                checked{++num1;}
            }
        }
        if (zLength <= 1)
            return;
        SyncLookupTableLength(zLength);
        int num3 = 0;
        while (num3 < yLength) {
            int num1 = 0;
            while (num1 < xLength) {
                int start = checked(num3 * inc2 + num1 * inc1);
                LinearFFT_Quick(data, start, inc3, zLength, direction);
                checked{++num1;}
            }
            checked{++num3;}
        }
    }
}


#endif //STRAIGHT_RAY_EXO_H
