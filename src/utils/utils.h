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
}
