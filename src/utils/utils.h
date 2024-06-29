namespace pptb::utils
{
    struct ascii_file_t
    {
        std::ifstream fh;
        std::string p_line{""};
        std::istringstream iss;
        std::string fn;
        std::size_t line_num = 0;
        bool f_eof = false;
    	ascii_file_t(const std::string& fname) : fh{fname}, fn{fname} {}
        bool eof() const { return f_eof; }
        
        std::size_t get_line_num() const { return line_num; }
        
        const std::string& next_line()
        {
            ++line_num;
            if (!std::getline(fh, p_line)) f_eof = true;
            iss.clear();
            iss.str(p_line);
            return p_line;
        }

        const std::string& line() const
        {
            return p_line;
        }
        
        template <typename data_t>
        void parse_sing(std::istringstream& iss_l, data_t& data)
        {
            iss_l >> data;
        }
        
        template <typename data_t>
        void r_parse(std::istringstream& iss_l, data_t& data)
        {
            parse_sing(iss_l, data);
        }
        
        void r_parse(std::istringstream& iss_l, std::string& data)
        {
            parse_sing(iss_l, data);
        }
        
        template <typename data_t, typename... datas_t>
		void r_parse(std::istringstream& iss_l, data_t& data, datas_t&... datas)
        {
            parse_sing(iss_l, data);
            r_parse(iss_l, datas...);
        }
        
        template <typename... datas_t>
        void parse(datas_t&... datas)
        {
            std::istringstream iss_loc(p_line);
            r_parse(iss_loc, datas...);
        }

        template <typename data_t>
        bool try_parse_sing(std::istringstream& iss_l, data_t& data)
        {
            iss_l >> data;
            return !iss_l.fail();
        }

        template <typename data_t>
        bool r_try_parse(std::istringstream& iss_l, data_t& data)
        {
            return try_parse_sing(iss_l, data);
        }
        
        template <typename data_t, typename... datas_t>
		bool r_try_parse(std::istringstream& iss_l, data_t& data, datas_t&... datas)
        {
            if (!try_parse_sing(iss_l, data)) return false;
            return r_try_parse(iss_l, datas...);
        }

        template <typename... datas_t>
        bool try_parse(datas_t&... datas)
        {
            return r_try_parse(iss, datas...);
        }
    };
}
