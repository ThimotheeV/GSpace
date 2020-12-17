#define CATCH_CONFIG_MAIN ArgParserTest
#include "catch.hpp"
#include <iostream>

#include "arg_parser.hpp"

TEST_CASE("arg_parser_file_test")
{
    SECTION("slice_and_lowercase")
    {
        REQUIRE(remove_spaces_tab_underscores("Text __with\tsome \t  whitespaces and_ _underscores") == "Textwithsomewhitespacesandunderscores");
        REQUIRE(remove_spaces_tab_underscores(str_tolower("Text __with\tsome \t  whitespaces and_ _underscores")) == "textwithsomewhitespacesandunderscores");
    }

    SECTION("slice_by_space")
    {
        auto result = slice_by_char(" 0201  003003 0102   0302 1011 01", ' ');

        REQUIRE(result.size() == 6);
        REQUIRE(result[0] == "0201");
        REQUIRE(result[1] == "003003");
        REQUIRE(result[2] == "0102");
        REQUIRE(result[3] == "0302");
        REQUIRE(result[4] == "1011");
        REQUIRE(result[5] == "01");
    }

    SECTION("slice_by_comma")
    {
        auto result = slice_by_char("Grange des Peres  ,  0201 003003 0102 0302 1011 01", ',');

        REQUIRE(result[0] == "Grange des Peres  ");
        REQUIRE(result[1] == "  0201 003003 0102 0302 1011 01");
    }

    SECTION("slice_by_equal")
    {
        auto result = slice_by_char("RecombinationRate=0.0000001", '=');

        REQUIRE(result[0] == "RecombinationRate");
        REQUIRE(result[1] == "0.0000001");
    }

    SECTION("slice_unix_file")
    {
        auto result = slice_unix_windows_file_by_line("line 1\nline 2");

        REQUIRE(result.size() == 2);
        REQUIRE(result[0] == "line 1");
        REQUIRE(result[1] == "line 2");
    }

    SECTION("slice_unix_file_with_comment")
    {
        auto result = slice_unix_windows_file_by_line("line 1\nline 2\n%line 3\nline 4\n%line 3");

        REQUIRE(result.size() == 3);
        REQUIRE(result[0] == "line 1");
        REQUIRE(result[1] == "line 2");
        REQUIRE(result[2] == "line 4");
    }

    SECTION("slice_windows_file")
    {
        auto result = slice_unix_windows_file_by_line("line 1\rline 2");

        REQUIRE(result.size() == 2);
        REQUIRE(result[0] == "line 1");
        REQUIRE(result[1] == "line 2");
    }

    SECTION("slice_by_string")
    {
        auto result = slice_str_vec_by_string({"mtDNA", "Pop", "Grange des Peres  ,  0201 003003 0102 0302 1011 01"}, "Pop");

        REQUIRE(result.size() == 2);
        REQUIRE(result[0][0] == "mtDNA");
        REQUIRE(result[1][0] == "Grange des Peres  ,  0201 003003 0102 0302 1011 01");
    }

    SECTION("key/value_splitting")
    {
        std::string txt{"truc=machin"};
        auto res = slice_by_char(txt, '=');
        REQUIRE(res[0] == "truc");
        REQUIRE(res[1] == "machin");

        std::string txt4{"truc=machin=bidule=35"};
        auto res2 = slice_by_char(txt4, '=');
        REQUIRE(res2[0] == "truc");
        REQUIRE(res2[1] == "machin");
        REQUIRE(res2[2] == "bidule");
        REQUIRE(res2[3] == "35");
    }
    SECTION("vector_parsing")
    {
        REQUIRE(int_vector_parse_by_comma("1,2,30") == std::vector<int>{1, 2, 30});
        REQUIRE(double_vector_parse_by_comma("10,-222,3.5") == std::vector<double>{10, -222, 3.5});
    }

    SECTION("bool convert_str_bool(std::string const &str)")
    {
        REQUIRE(convert_str_bool("true"));
        REQUIRE(convert_str_bool("t"));
        REQUIRE(!convert_str_bool("false"));
        REQUIRE(!convert_str_bool("f"));
        REQUIRE_THROWS(convert_str_bool("Tru"));
        REQUIRE_THROWS(convert_str_bool("Falses"));
    }
}
