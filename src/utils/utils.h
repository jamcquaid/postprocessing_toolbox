namespace pptb::utils
{
	struct file_directory_t
	{
		std::size_t nFiles;
		std::vector<std::string> filenames;

		file_directory_t()
		{
			print("Parsing provided directory...");

			// Build command string
			auto cmd_out = std::system("ls data > files.log");

			// Open file
			spade::utils::ascii_file_t file("files.log");

			// Loop until EOF
			nFiles = 0;
			while (!file.eof())
			{
				// Next line
				file.next_line();

				// EOF?
				if (file.eof()) break;

				// Store filename
				filenames.push_back("data/"+file.p_line);

				// Count files
				nFiles++;
			}

			print("Number of files found =  ",nFiles);
		}
	};

	std::vector<std::string> splitString(const std::string& str, const std::string& delimiter)
	{
		std::vector<std::string> tokens;
		size_t start = 0;
		size_t end = str.find(delimiter);
		while (end != std::string::npos) {
			tokens.push_back(str.substr(start, end - start));
			start = end + delimiter.length();
			end = str.find(delimiter, start);
		}
		tokens.push_back(str.substr(start));
		return tokens;
	}
}
