#pragma once

// #include <algorithm>
#include <vector>
#include <cstdint>
#include <cstring>
#include <iostream>
#include <map> 
#include <mutex>
#include <condition_variable> 
using namespace std;

class File_Handle_2 {
private:
	bool input_mode;
	FILE* f;
	size_t f_offset;
	string file_name;

	struct part_t{
		size_t offset;
		size_t size;

		part_t() : offset(0), size(0)
		{};

		part_t(size_t _offset, size_t _size) : offset(_offset), size(_size)
		{};
	};

	typedef struct {
		string stream_name;
		size_t cur_id;
		vector<part_t> parts;
	} stream_t;

	map<int, stream_t> m_streams;
	mutex mtx;

	bool serialize();
	bool deserialize();
	size_t WriteFixed(size_t x, FILE* file);
	size_t Write(size_t x, FILE* file);
	size_t Write(string s, FILE* file);
	size_t read_fixed(size_t& x, FILE* file);
	size_t read(size_t& x, FILE* file);
	size_t read(string& s, FILE* file);

public:
	File_Handle_2(bool _input_mode);
	~File_Handle_2();

	bool Open(string _file_name);
	bool Close();

	int RegisterStream(string stream_name);
	int GetStreamId(string stream_name);

	bool AddPart(int stream_id, vector<uint8_t> &v_data, size_t metadata = 0);
	int AddPartPrepare(int stream_id);
	bool AddPartComplete(int stream_id, int part_id, vector<uint8_t>& v_data);
	void AddParamsPart(int stream_id,vector<uint8_t> & v_data);

	bool GetPart(int stream_id, vector<uint8_t> &v_data);
	void SetRawSize(int stream_id, size_t raw_size);
	size_t GetRawSize(int stream_id);
	size_t GetCompressedSize(int stream_id);
	bool ResetStreamPartIterator(int stream_id);

	bool LinkStream(int stream_id, string stream_name, int target_id);


	size_t GetNoStreams()
	{
		lock_guard<mutex> lck(mtx);

		return m_streams.size();
	}
};