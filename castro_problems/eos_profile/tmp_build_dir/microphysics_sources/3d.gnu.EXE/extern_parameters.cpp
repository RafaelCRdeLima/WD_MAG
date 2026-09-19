#include <extern_parameters.H>
#include <AMReX_ParmParse.H>

#include <AMReX_REAL.H>

  namespace eos_rp {
    AMREX_GPU_MANAGED bool use_eos_coulomb;
    AMREX_GPU_MANAGED bool eos_input_is_constant;
    AMREX_GPU_MANAGED amrex::Real eos_ttol;
    AMREX_GPU_MANAGED amrex::Real eos_dtol;
    AMREX_GPU_MANAGED amrex::Real prad_limiter_rho_c;
    AMREX_GPU_MANAGED amrex::Real prad_limiter_delta_rho;
  }
  namespace network_rp {
    AMREX_GPU_MANAGED amrex::Real small_x;
    AMREX_GPU_MANAGED bool use_tables;
    AMREX_GPU_MANAGED bool use_c12ag_deboer17;
  }
  namespace unit_test_rp {
    std::string primary_species_1;
    std::string primary_species_2;
    std::string primary_species_3;
    AMREX_GPU_MANAGED amrex::Real X1;
    AMREX_GPU_MANAGED amrex::Real X2;
    AMREX_GPU_MANAGED amrex::Real X3;
    AMREX_GPU_MANAGED amrex::Real X4;
    AMREX_GPU_MANAGED amrex::Real X5;
    AMREX_GPU_MANAGED amrex::Real X6;
    AMREX_GPU_MANAGED amrex::Real X7;
    AMREX_GPU_MANAGED amrex::Real X8;
    AMREX_GPU_MANAGED amrex::Real X9;
    AMREX_GPU_MANAGED amrex::Real X10;
    AMREX_GPU_MANAGED amrex::Real X11;
    AMREX_GPU_MANAGED amrex::Real X12;
    AMREX_GPU_MANAGED amrex::Real X13;
    AMREX_GPU_MANAGED amrex::Real X14;
    AMREX_GPU_MANAGED amrex::Real X15;
    AMREX_GPU_MANAGED amrex::Real X16;
    AMREX_GPU_MANAGED amrex::Real X17;
    AMREX_GPU_MANAGED amrex::Real X18;
    AMREX_GPU_MANAGED amrex::Real X19;
    AMREX_GPU_MANAGED amrex::Real X20;
    AMREX_GPU_MANAGED amrex::Real X21;
    AMREX_GPU_MANAGED amrex::Real X22;
    AMREX_GPU_MANAGED amrex::Real X23;
    AMREX_GPU_MANAGED amrex::Real X24;
    AMREX_GPU_MANAGED amrex::Real X25;
    AMREX_GPU_MANAGED amrex::Real X26;
    AMREX_GPU_MANAGED amrex::Real X27;
    AMREX_GPU_MANAGED amrex::Real X28;
    AMREX_GPU_MANAGED amrex::Real X29;
    AMREX_GPU_MANAGED amrex::Real X30;
    AMREX_GPU_MANAGED amrex::Real X31;
    AMREX_GPU_MANAGED amrex::Real X32;
    AMREX_GPU_MANAGED amrex::Real X33;
    AMREX_GPU_MANAGED amrex::Real X34;
    AMREX_GPU_MANAGED amrex::Real X35;
    AMREX_GPU_MANAGED amrex::Real X36;
    AMREX_GPU_MANAGED amrex::Real X37;
    AMREX_GPU_MANAGED amrex::Real X38;
    AMREX_GPU_MANAGED amrex::Real X39;
    AMREX_GPU_MANAGED amrex::Real X40;
    AMREX_GPU_MANAGED amrex::Real X41;
    AMREX_GPU_MANAGED amrex::Real X42;
    AMREX_GPU_MANAGED amrex::Real X43;
    AMREX_GPU_MANAGED amrex::Real X44;
    AMREX_GPU_MANAGED amrex::Real X45;
    AMREX_GPU_MANAGED amrex::Real X46;
    AMREX_GPU_MANAGED amrex::Real X47;
    AMREX_GPU_MANAGED amrex::Real X48;
    AMREX_GPU_MANAGED amrex::Real X49;
    AMREX_GPU_MANAGED amrex::Real X50;
    AMREX_GPU_MANAGED amrex::Real X51;
    AMREX_GPU_MANAGED amrex::Real X52;
    AMREX_GPU_MANAGED amrex::Real X53;
    AMREX_GPU_MANAGED amrex::Real X54;
    AMREX_GPU_MANAGED amrex::Real X55;
    AMREX_GPU_MANAGED amrex::Real X56;
    AMREX_GPU_MANAGED amrex::Real X57;
    AMREX_GPU_MANAGED amrex::Real X58;
    AMREX_GPU_MANAGED amrex::Real X59;
    AMREX_GPU_MANAGED amrex::Real X60;
    AMREX_GPU_MANAGED amrex::Real X61;
    AMREX_GPU_MANAGED amrex::Real X62;
    AMREX_GPU_MANAGED amrex::Real X63;
    AMREX_GPU_MANAGED amrex::Real X64;
    AMREX_GPU_MANAGED amrex::Real X65;
    AMREX_GPU_MANAGED amrex::Real X66;
    AMREX_GPU_MANAGED amrex::Real X67;
    AMREX_GPU_MANAGED amrex::Real X68;
    AMREX_GPU_MANAGED amrex::Real X69;
    AMREX_GPU_MANAGED amrex::Real X70;
    AMREX_GPU_MANAGED amrex::Real X71;
    AMREX_GPU_MANAGED amrex::Real X72;
    AMREX_GPU_MANAGED amrex::Real X73;
    AMREX_GPU_MANAGED amrex::Real X74;
    AMREX_GPU_MANAGED amrex::Real X75;
    AMREX_GPU_MANAGED amrex::Real X76;
    AMREX_GPU_MANAGED amrex::Real X77;
    AMREX_GPU_MANAGED amrex::Real X78;
    AMREX_GPU_MANAGED amrex::Real X79;
    AMREX_GPU_MANAGED amrex::Real X80;
    AMREX_GPU_MANAGED amrex::Real X81;
    AMREX_GPU_MANAGED amrex::Real X82;
    AMREX_GPU_MANAGED amrex::Real X83;
    AMREX_GPU_MANAGED amrex::Real X84;
    AMREX_GPU_MANAGED amrex::Real X85;
    AMREX_GPU_MANAGED amrex::Real X86;
    AMREX_GPU_MANAGED amrex::Real X87;
    AMREX_GPU_MANAGED amrex::Real X88;
    AMREX_GPU_MANAGED amrex::Real X89;
    AMREX_GPU_MANAGED amrex::Real X90;
    AMREX_GPU_MANAGED amrex::Real X91;
    AMREX_GPU_MANAGED amrex::Real X92;
    AMREX_GPU_MANAGED amrex::Real X93;
    AMREX_GPU_MANAGED amrex::Real X94;
    AMREX_GPU_MANAGED amrex::Real X95;
    AMREX_GPU_MANAGED amrex::Real X96;
    AMREX_GPU_MANAGED amrex::Real X97;
    AMREX_GPU_MANAGED amrex::Real X98;
    AMREX_GPU_MANAGED amrex::Real X99;
    AMREX_GPU_MANAGED amrex::Real X100;
    AMREX_GPU_MANAGED amrex::Real X101;
    AMREX_GPU_MANAGED amrex::Real X102;
    AMREX_GPU_MANAGED amrex::Real X103;
    AMREX_GPU_MANAGED amrex::Real X104;
    AMREX_GPU_MANAGED amrex::Real X105;
    AMREX_GPU_MANAGED amrex::Real X106;
    AMREX_GPU_MANAGED amrex::Real X107;
    AMREX_GPU_MANAGED amrex::Real X108;
    AMREX_GPU_MANAGED amrex::Real X109;
    AMREX_GPU_MANAGED amrex::Real X110;
    AMREX_GPU_MANAGED amrex::Real X111;
    AMREX_GPU_MANAGED amrex::Real X112;
    AMREX_GPU_MANAGED amrex::Real X113;
    AMREX_GPU_MANAGED amrex::Real X114;
    AMREX_GPU_MANAGED amrex::Real X115;
    AMREX_GPU_MANAGED amrex::Real X116;
    AMREX_GPU_MANAGED amrex::Real X117;
    AMREX_GPU_MANAGED amrex::Real X118;
    AMREX_GPU_MANAGED amrex::Real X119;
    AMREX_GPU_MANAGED amrex::Real X120;
    AMREX_GPU_MANAGED amrex::Real X121;
    AMREX_GPU_MANAGED amrex::Real X122;
    AMREX_GPU_MANAGED amrex::Real X123;
    AMREX_GPU_MANAGED amrex::Real X124;
    AMREX_GPU_MANAGED amrex::Real X125;
    AMREX_GPU_MANAGED amrex::Real X126;
    AMREX_GPU_MANAGED amrex::Real X127;
    AMREX_GPU_MANAGED amrex::Real X128;
    AMREX_GPU_MANAGED amrex::Real X129;
    AMREX_GPU_MANAGED amrex::Real X130;
    AMREX_GPU_MANAGED amrex::Real X131;
    AMREX_GPU_MANAGED amrex::Real X132;
    AMREX_GPU_MANAGED amrex::Real X133;
    AMREX_GPU_MANAGED amrex::Real X134;
    AMREX_GPU_MANAGED amrex::Real X135;
    AMREX_GPU_MANAGED amrex::Real X136;
    AMREX_GPU_MANAGED amrex::Real X137;
    AMREX_GPU_MANAGED amrex::Real X138;
    AMREX_GPU_MANAGED amrex::Real X139;
    AMREX_GPU_MANAGED amrex::Real X140;
    AMREX_GPU_MANAGED amrex::Real X141;
    AMREX_GPU_MANAGED amrex::Real X142;
    AMREX_GPU_MANAGED amrex::Real X143;
    AMREX_GPU_MANAGED amrex::Real X144;
    AMREX_GPU_MANAGED amrex::Real X145;
    AMREX_GPU_MANAGED amrex::Real X146;
    AMREX_GPU_MANAGED amrex::Real X147;
    AMREX_GPU_MANAGED amrex::Real X148;
    AMREX_GPU_MANAGED amrex::Real X149;
    AMREX_GPU_MANAGED amrex::Real X150;
    AMREX_GPU_MANAGED amrex::Real X151;
    AMREX_GPU_MANAGED amrex::Real X152;
    AMREX_GPU_MANAGED amrex::Real X153;
    AMREX_GPU_MANAGED amrex::Real X154;
    AMREX_GPU_MANAGED amrex::Real X155;
    AMREX_GPU_MANAGED amrex::Real X156;
    AMREX_GPU_MANAGED amrex::Real X157;
    AMREX_GPU_MANAGED amrex::Real X158;
    AMREX_GPU_MANAGED amrex::Real X159;
    AMREX_GPU_MANAGED amrex::Real X160;
    AMREX_GPU_MANAGED amrex::Real X161;
    AMREX_GPU_MANAGED amrex::Real X162;
    AMREX_GPU_MANAGED amrex::Real X163;
    AMREX_GPU_MANAGED amrex::Real X164;
    AMREX_GPU_MANAGED amrex::Real X165;
    AMREX_GPU_MANAGED amrex::Real X166;
    AMREX_GPU_MANAGED amrex::Real X167;
    AMREX_GPU_MANAGED amrex::Real X168;
    AMREX_GPU_MANAGED amrex::Real X169;
    AMREX_GPU_MANAGED amrex::Real X170;
    AMREX_GPU_MANAGED amrex::Real X171;
    AMREX_GPU_MANAGED amrex::Real X172;
    AMREX_GPU_MANAGED amrex::Real X173;
    AMREX_GPU_MANAGED amrex::Real X174;
    AMREX_GPU_MANAGED amrex::Real X175;
    AMREX_GPU_MANAGED amrex::Real X176;
    AMREX_GPU_MANAGED amrex::Real X177;
    AMREX_GPU_MANAGED amrex::Real X178;
    AMREX_GPU_MANAGED amrex::Real X179;
    AMREX_GPU_MANAGED amrex::Real X180;
    AMREX_GPU_MANAGED amrex::Real X181;
    AMREX_GPU_MANAGED amrex::Real X182;
    AMREX_GPU_MANAGED amrex::Real X183;
    AMREX_GPU_MANAGED amrex::Real X184;
    AMREX_GPU_MANAGED amrex::Real X185;
    AMREX_GPU_MANAGED amrex::Real X186;
    AMREX_GPU_MANAGED amrex::Real X187;
    AMREX_GPU_MANAGED amrex::Real X188;
    AMREX_GPU_MANAGED amrex::Real X189;
    AMREX_GPU_MANAGED amrex::Real X190;
    AMREX_GPU_MANAGED amrex::Real X191;
    AMREX_GPU_MANAGED amrex::Real X192;
    AMREX_GPU_MANAGED amrex::Real X193;
    AMREX_GPU_MANAGED amrex::Real X194;
    AMREX_GPU_MANAGED amrex::Real X195;
    AMREX_GPU_MANAGED amrex::Real X196;
    AMREX_GPU_MANAGED amrex::Real X197;
    AMREX_GPU_MANAGED amrex::Real X198;
    AMREX_GPU_MANAGED amrex::Real X199;
    AMREX_GPU_MANAGED amrex::Real X200;
    AMREX_GPU_MANAGED amrex::Real X201;
    AMREX_GPU_MANAGED amrex::Real X202;
    AMREX_GPU_MANAGED amrex::Real X203;
    AMREX_GPU_MANAGED amrex::Real X204;
    AMREX_GPU_MANAGED amrex::Real X205;
    AMREX_GPU_MANAGED amrex::Real X206;
    AMREX_GPU_MANAGED amrex::Real X207;
    AMREX_GPU_MANAGED amrex::Real X208;
    AMREX_GPU_MANAGED amrex::Real X209;
    AMREX_GPU_MANAGED amrex::Real X210;
    AMREX_GPU_MANAGED amrex::Real X211;
    AMREX_GPU_MANAGED amrex::Real X212;
    AMREX_GPU_MANAGED amrex::Real X213;
    AMREX_GPU_MANAGED amrex::Real X214;
    AMREX_GPU_MANAGED amrex::Real X215;
    AMREX_GPU_MANAGED amrex::Real X216;
    AMREX_GPU_MANAGED amrex::Real X217;
    AMREX_GPU_MANAGED amrex::Real X218;
    AMREX_GPU_MANAGED amrex::Real X219;
    AMREX_GPU_MANAGED amrex::Real X220;
    AMREX_GPU_MANAGED amrex::Real X221;
    AMREX_GPU_MANAGED amrex::Real X222;
    AMREX_GPU_MANAGED amrex::Real X223;
    AMREX_GPU_MANAGED amrex::Real X224;
    AMREX_GPU_MANAGED amrex::Real X225;
    AMREX_GPU_MANAGED amrex::Real X226;
    AMREX_GPU_MANAGED amrex::Real X227;
    AMREX_GPU_MANAGED amrex::Real X228;
    AMREX_GPU_MANAGED amrex::Real X229;
    AMREX_GPU_MANAGED amrex::Real X230;
    AMREX_GPU_MANAGED amrex::Real X231;
    AMREX_GPU_MANAGED amrex::Real X232;
    AMREX_GPU_MANAGED amrex::Real X233;
    AMREX_GPU_MANAGED amrex::Real X234;
    AMREX_GPU_MANAGED amrex::Real X235;
    AMREX_GPU_MANAGED amrex::Real X236;
    AMREX_GPU_MANAGED amrex::Real X237;
    AMREX_GPU_MANAGED amrex::Real X238;
    AMREX_GPU_MANAGED amrex::Real X239;
    AMREX_GPU_MANAGED amrex::Real X240;
    AMREX_GPU_MANAGED amrex::Real X241;
    AMREX_GPU_MANAGED amrex::Real X242;
    AMREX_GPU_MANAGED amrex::Real X243;
    AMREX_GPU_MANAGED amrex::Real X244;
    AMREX_GPU_MANAGED amrex::Real X245;
    AMREX_GPU_MANAGED amrex::Real X246;
    AMREX_GPU_MANAGED amrex::Real X247;
    AMREX_GPU_MANAGED amrex::Real X248;
    AMREX_GPU_MANAGED amrex::Real X249;
    AMREX_GPU_MANAGED amrex::Real X250;
    AMREX_GPU_MANAGED amrex::Real X251;
    AMREX_GPU_MANAGED amrex::Real X252;
    AMREX_GPU_MANAGED amrex::Real X253;
    AMREX_GPU_MANAGED amrex::Real X254;
    AMREX_GPU_MANAGED amrex::Real X255;
    AMREX_GPU_MANAGED amrex::Real X256;
    AMREX_GPU_MANAGED amrex::Real X257;
    AMREX_GPU_MANAGED amrex::Real X258;
    AMREX_GPU_MANAGED amrex::Real X259;
    AMREX_GPU_MANAGED amrex::Real X260;
    AMREX_GPU_MANAGED amrex::Real X261;
    AMREX_GPU_MANAGED amrex::Real X262;
    AMREX_GPU_MANAGED amrex::Real X263;
    AMREX_GPU_MANAGED amrex::Real X264;
    AMREX_GPU_MANAGED amrex::Real X265;
    AMREX_GPU_MANAGED amrex::Real X266;
    AMREX_GPU_MANAGED amrex::Real X267;
    AMREX_GPU_MANAGED amrex::Real X268;
    AMREX_GPU_MANAGED amrex::Real X269;
    AMREX_GPU_MANAGED amrex::Real X270;
    AMREX_GPU_MANAGED amrex::Real X271;
    AMREX_GPU_MANAGED amrex::Real X272;
    AMREX_GPU_MANAGED amrex::Real X273;
    AMREX_GPU_MANAGED amrex::Real X274;
    AMREX_GPU_MANAGED amrex::Real X275;
    AMREX_GPU_MANAGED amrex::Real X276;
    AMREX_GPU_MANAGED amrex::Real X277;
    AMREX_GPU_MANAGED amrex::Real X278;
    AMREX_GPU_MANAGED amrex::Real X279;
    AMREX_GPU_MANAGED amrex::Real X280;
    AMREX_GPU_MANAGED amrex::Real X281;
    AMREX_GPU_MANAGED amrex::Real X282;
    AMREX_GPU_MANAGED amrex::Real X283;
    AMREX_GPU_MANAGED amrex::Real X284;
    AMREX_GPU_MANAGED amrex::Real X285;
    AMREX_GPU_MANAGED amrex::Real X286;
    AMREX_GPU_MANAGED amrex::Real X287;
    AMREX_GPU_MANAGED amrex::Real X288;
    AMREX_GPU_MANAGED amrex::Real X289;
    AMREX_GPU_MANAGED amrex::Real X290;
    AMREX_GPU_MANAGED amrex::Real X291;
    AMREX_GPU_MANAGED amrex::Real X292;
    AMREX_GPU_MANAGED amrex::Real X293;
    AMREX_GPU_MANAGED amrex::Real X294;
    AMREX_GPU_MANAGED amrex::Real X295;
    AMREX_GPU_MANAGED amrex::Real X296;
    AMREX_GPU_MANAGED amrex::Real X297;
    AMREX_GPU_MANAGED amrex::Real X298;
    AMREX_GPU_MANAGED amrex::Real X299;
    AMREX_GPU_MANAGED amrex::Real X300;
    AMREX_GPU_MANAGED amrex::Real X301;
    AMREX_GPU_MANAGED amrex::Real X302;
    AMREX_GPU_MANAGED amrex::Real X303;
    AMREX_GPU_MANAGED amrex::Real X304;
    AMREX_GPU_MANAGED amrex::Real X305;
    AMREX_GPU_MANAGED amrex::Real X306;
    AMREX_GPU_MANAGED amrex::Real X307;
    AMREX_GPU_MANAGED amrex::Real X308;
    AMREX_GPU_MANAGED amrex::Real X309;
    AMREX_GPU_MANAGED amrex::Real X310;
    AMREX_GPU_MANAGED amrex::Real X311;
    AMREX_GPU_MANAGED amrex::Real X312;
    AMREX_GPU_MANAGED amrex::Real X313;
    AMREX_GPU_MANAGED amrex::Real X314;
    AMREX_GPU_MANAGED amrex::Real X315;
    AMREX_GPU_MANAGED amrex::Real X316;
    AMREX_GPU_MANAGED amrex::Real X317;
    AMREX_GPU_MANAGED amrex::Real X318;
    AMREX_GPU_MANAGED amrex::Real X319;
    AMREX_GPU_MANAGED amrex::Real X320;
    AMREX_GPU_MANAGED amrex::Real X321;
    AMREX_GPU_MANAGED amrex::Real X322;
    AMREX_GPU_MANAGED amrex::Real X323;
    AMREX_GPU_MANAGED amrex::Real X324;
    AMREX_GPU_MANAGED amrex::Real X325;
    AMREX_GPU_MANAGED amrex::Real X326;
    AMREX_GPU_MANAGED amrex::Real X327;
    AMREX_GPU_MANAGED amrex::Real X328;
    AMREX_GPU_MANAGED amrex::Real X329;
    AMREX_GPU_MANAGED amrex::Real X330;
    AMREX_GPU_MANAGED amrex::Real X331;
    AMREX_GPU_MANAGED amrex::Real X332;
    AMREX_GPU_MANAGED amrex::Real X333;
    AMREX_GPU_MANAGED amrex::Real X334;
    AMREX_GPU_MANAGED amrex::Real X335;
    AMREX_GPU_MANAGED amrex::Real X336;
    AMREX_GPU_MANAGED amrex::Real X337;
    AMREX_GPU_MANAGED amrex::Real X338;
    AMREX_GPU_MANAGED amrex::Real X339;
    AMREX_GPU_MANAGED amrex::Real X340;
    AMREX_GPU_MANAGED amrex::Real X341;
    AMREX_GPU_MANAGED amrex::Real X342;
    AMREX_GPU_MANAGED amrex::Real X343;
    AMREX_GPU_MANAGED amrex::Real X344;
    AMREX_GPU_MANAGED amrex::Real X345;
    AMREX_GPU_MANAGED amrex::Real X346;
    AMREX_GPU_MANAGED amrex::Real X347;
    AMREX_GPU_MANAGED amrex::Real X348;
    AMREX_GPU_MANAGED amrex::Real X349;
    AMREX_GPU_MANAGED amrex::Real X350;
    AMREX_GPU_MANAGED amrex::Real X351;
    AMREX_GPU_MANAGED amrex::Real X352;
    AMREX_GPU_MANAGED amrex::Real X353;
    AMREX_GPU_MANAGED amrex::Real X354;
    AMREX_GPU_MANAGED amrex::Real X355;
    AMREX_GPU_MANAGED amrex::Real X356;
    AMREX_GPU_MANAGED amrex::Real X357;
    AMREX_GPU_MANAGED amrex::Real X358;
    AMREX_GPU_MANAGED amrex::Real X359;
    AMREX_GPU_MANAGED amrex::Real X360;
    AMREX_GPU_MANAGED amrex::Real X361;
    AMREX_GPU_MANAGED amrex::Real X362;
    AMREX_GPU_MANAGED amrex::Real X363;
    AMREX_GPU_MANAGED amrex::Real X364;
    AMREX_GPU_MANAGED amrex::Real X365;
    AMREX_GPU_MANAGED amrex::Real X366;
    AMREX_GPU_MANAGED amrex::Real X367;
    AMREX_GPU_MANAGED amrex::Real X368;
    AMREX_GPU_MANAGED amrex::Real X369;
    AMREX_GPU_MANAGED amrex::Real X370;
    AMREX_GPU_MANAGED amrex::Real X371;
    AMREX_GPU_MANAGED amrex::Real X372;
    AMREX_GPU_MANAGED amrex::Real X373;
    AMREX_GPU_MANAGED amrex::Real X374;
    AMREX_GPU_MANAGED amrex::Real X375;
    AMREX_GPU_MANAGED amrex::Real X376;
    AMREX_GPU_MANAGED amrex::Real X377;
    AMREX_GPU_MANAGED amrex::Real X378;
    AMREX_GPU_MANAGED amrex::Real X379;
    AMREX_GPU_MANAGED amrex::Real X380;
    AMREX_GPU_MANAGED amrex::Real X381;
    AMREX_GPU_MANAGED amrex::Real X382;
    AMREX_GPU_MANAGED amrex::Real X383;
    AMREX_GPU_MANAGED amrex::Real X384;
    AMREX_GPU_MANAGED amrex::Real X385;
    AMREX_GPU_MANAGED amrex::Real X386;
    AMREX_GPU_MANAGED amrex::Real X387;
    AMREX_GPU_MANAGED amrex::Real X388;
    AMREX_GPU_MANAGED amrex::Real X389;
    AMREX_GPU_MANAGED amrex::Real X390;
    AMREX_GPU_MANAGED amrex::Real X391;
    AMREX_GPU_MANAGED amrex::Real X392;
    AMREX_GPU_MANAGED amrex::Real X393;
    AMREX_GPU_MANAGED amrex::Real X394;
    AMREX_GPU_MANAGED amrex::Real X395;
    AMREX_GPU_MANAGED amrex::Real X396;
    AMREX_GPU_MANAGED amrex::Real X397;
    AMREX_GPU_MANAGED amrex::Real X398;
    AMREX_GPU_MANAGED amrex::Real X399;
    AMREX_GPU_MANAGED amrex::Real X400;
    AMREX_GPU_MANAGED bool uniform_xn;
    AMREX_GPU_MANAGED amrex::Real small_temp;
    AMREX_GPU_MANAGED amrex::Real small_dens;
  }

  extern_t init_extern_parameters() {
    using namespace amrex::literals;

    extern_t params;

    // get the value from the inputs file
    {
      const amrex::ParmParse pp("eos");
      eos_rp::use_eos_coulomb = true;
      pp.query("use_eos_coulomb", eos_rp::use_eos_coulomb);
      pp.query("use_eos_coulomb", params.eos.use_eos_coulomb);

      eos_rp::eos_input_is_constant = true;
      pp.query("eos_input_is_constant", eos_rp::eos_input_is_constant);
      pp.query("eos_input_is_constant", params.eos.eos_input_is_constant);

      eos_rp::eos_ttol = 1.0e-8_rt;
      pp.query("eos_ttol", eos_rp::eos_ttol);
      pp.query("eos_ttol", params.eos.eos_ttol);

      eos_rp::eos_dtol = 1.0e-8_rt;
      pp.query("eos_dtol", eos_rp::eos_dtol);
      pp.query("eos_dtol", params.eos.eos_dtol);

      eos_rp::prad_limiter_rho_c = -1.0e0_rt;
      pp.query("prad_limiter_rho_c", eos_rp::prad_limiter_rho_c);
      pp.query("prad_limiter_rho_c", params.eos.prad_limiter_rho_c);

      eos_rp::prad_limiter_delta_rho = -1.0e0_rt;
      pp.query("prad_limiter_delta_rho", eos_rp::prad_limiter_delta_rho);
      pp.query("prad_limiter_delta_rho", params.eos.prad_limiter_delta_rho);

    }
    {
      const amrex::ParmParse pp("network");
      network_rp::small_x = 1.e-30_rt;
      pp.query("small_x", network_rp::small_x);
      pp.query("small_x", params.network.small_x);

      network_rp::use_tables = false;
      pp.query("use_tables", network_rp::use_tables);
      pp.query("use_tables", params.network.use_tables);

      network_rp::use_c12ag_deboer17 = false;
      pp.query("use_c12ag_deboer17", network_rp::use_c12ag_deboer17);
      pp.query("use_c12ag_deboer17", params.network.use_c12ag_deboer17);

    }
    {
      const amrex::ParmParse pp("unit_test");
      unit_test_rp::primary_species_1 = "";
      pp.query("primary_species_1", unit_test_rp::primary_species_1);
      pp.query("primary_species_1", params.unit_test.primary_species_1);

      unit_test_rp::primary_species_2 = "";
      pp.query("primary_species_2", unit_test_rp::primary_species_2);
      pp.query("primary_species_2", params.unit_test.primary_species_2);

      unit_test_rp::primary_species_3 = "";
      pp.query("primary_species_3", unit_test_rp::primary_species_3);
      pp.query("primary_species_3", params.unit_test.primary_species_3);

      unit_test_rp::X1 = 1.0e0_rt;
      pp.query("X1", unit_test_rp::X1);
      pp.query("X1", params.unit_test.X1);

      unit_test_rp::X2 = 0.0e0_rt;
      pp.query("X2", unit_test_rp::X2);
      pp.query("X2", params.unit_test.X2);

      unit_test_rp::X3 = 0.0e0_rt;
      pp.query("X3", unit_test_rp::X3);
      pp.query("X3", params.unit_test.X3);

      unit_test_rp::X4 = 0.0e0_rt;
      pp.query("X4", unit_test_rp::X4);
      pp.query("X4", params.unit_test.X4);

      unit_test_rp::X5 = 0.0e0_rt;
      pp.query("X5", unit_test_rp::X5);
      pp.query("X5", params.unit_test.X5);

      unit_test_rp::X6 = 0.0e0_rt;
      pp.query("X6", unit_test_rp::X6);
      pp.query("X6", params.unit_test.X6);

      unit_test_rp::X7 = 0.0e0_rt;
      pp.query("X7", unit_test_rp::X7);
      pp.query("X7", params.unit_test.X7);

      unit_test_rp::X8 = 0.0e0_rt;
      pp.query("X8", unit_test_rp::X8);
      pp.query("X8", params.unit_test.X8);

      unit_test_rp::X9 = 0.0e0_rt;
      pp.query("X9", unit_test_rp::X9);
      pp.query("X9", params.unit_test.X9);

      unit_test_rp::X10 = 0.0e0_rt;
      pp.query("X10", unit_test_rp::X10);
      pp.query("X10", params.unit_test.X10);

      unit_test_rp::X11 = 0.0e0_rt;
      pp.query("X11", unit_test_rp::X11);
      pp.query("X11", params.unit_test.X11);

      unit_test_rp::X12 = 0.0e0_rt;
      pp.query("X12", unit_test_rp::X12);
      pp.query("X12", params.unit_test.X12);

      unit_test_rp::X13 = 0.0e0_rt;
      pp.query("X13", unit_test_rp::X13);
      pp.query("X13", params.unit_test.X13);

      unit_test_rp::X14 = 0.0e0_rt;
      pp.query("X14", unit_test_rp::X14);
      pp.query("X14", params.unit_test.X14);

      unit_test_rp::X15 = 0.0e0_rt;
      pp.query("X15", unit_test_rp::X15);
      pp.query("X15", params.unit_test.X15);

      unit_test_rp::X16 = 0.0e0_rt;
      pp.query("X16", unit_test_rp::X16);
      pp.query("X16", params.unit_test.X16);

      unit_test_rp::X17 = 0.0e0_rt;
      pp.query("X17", unit_test_rp::X17);
      pp.query("X17", params.unit_test.X17);

      unit_test_rp::X18 = 0.0e0_rt;
      pp.query("X18", unit_test_rp::X18);
      pp.query("X18", params.unit_test.X18);

      unit_test_rp::X19 = 0.0e0_rt;
      pp.query("X19", unit_test_rp::X19);
      pp.query("X19", params.unit_test.X19);

      unit_test_rp::X20 = 0.0e0_rt;
      pp.query("X20", unit_test_rp::X20);
      pp.query("X20", params.unit_test.X20);

      unit_test_rp::X21 = 0.0e0_rt;
      pp.query("X21", unit_test_rp::X21);
      pp.query("X21", params.unit_test.X21);

      unit_test_rp::X22 = 0.0e0_rt;
      pp.query("X22", unit_test_rp::X22);
      pp.query("X22", params.unit_test.X22);

      unit_test_rp::X23 = 0.0e0_rt;
      pp.query("X23", unit_test_rp::X23);
      pp.query("X23", params.unit_test.X23);

      unit_test_rp::X24 = 0.0e0_rt;
      pp.query("X24", unit_test_rp::X24);
      pp.query("X24", params.unit_test.X24);

      unit_test_rp::X25 = 0.0e0_rt;
      pp.query("X25", unit_test_rp::X25);
      pp.query("X25", params.unit_test.X25);

      unit_test_rp::X26 = 0.0e0_rt;
      pp.query("X26", unit_test_rp::X26);
      pp.query("X26", params.unit_test.X26);

      unit_test_rp::X27 = 0.0e0_rt;
      pp.query("X27", unit_test_rp::X27);
      pp.query("X27", params.unit_test.X27);

      unit_test_rp::X28 = 0.0e0_rt;
      pp.query("X28", unit_test_rp::X28);
      pp.query("X28", params.unit_test.X28);

      unit_test_rp::X29 = 0.0e0_rt;
      pp.query("X29", unit_test_rp::X29);
      pp.query("X29", params.unit_test.X29);

      unit_test_rp::X30 = 0.0e0_rt;
      pp.query("X30", unit_test_rp::X30);
      pp.query("X30", params.unit_test.X30);

      unit_test_rp::X31 = 0.0e0_rt;
      pp.query("X31", unit_test_rp::X31);
      pp.query("X31", params.unit_test.X31);

      unit_test_rp::X32 = 0.0e0_rt;
      pp.query("X32", unit_test_rp::X32);
      pp.query("X32", params.unit_test.X32);

      unit_test_rp::X33 = 0.0e0_rt;
      pp.query("X33", unit_test_rp::X33);
      pp.query("X33", params.unit_test.X33);

      unit_test_rp::X34 = 0.0e0_rt;
      pp.query("X34", unit_test_rp::X34);
      pp.query("X34", params.unit_test.X34);

      unit_test_rp::X35 = 0.0e0_rt;
      pp.query("X35", unit_test_rp::X35);
      pp.query("X35", params.unit_test.X35);

      unit_test_rp::X36 = 0.0e0_rt;
      pp.query("X36", unit_test_rp::X36);
      pp.query("X36", params.unit_test.X36);

      unit_test_rp::X37 = 0.0e0_rt;
      pp.query("X37", unit_test_rp::X37);
      pp.query("X37", params.unit_test.X37);

      unit_test_rp::X38 = 0.0e0_rt;
      pp.query("X38", unit_test_rp::X38);
      pp.query("X38", params.unit_test.X38);

      unit_test_rp::X39 = 0.0e0_rt;
      pp.query("X39", unit_test_rp::X39);
      pp.query("X39", params.unit_test.X39);

      unit_test_rp::X40 = 0.0e0_rt;
      pp.query("X40", unit_test_rp::X40);
      pp.query("X40", params.unit_test.X40);

      unit_test_rp::X41 = 0.0e0_rt;
      pp.query("X41", unit_test_rp::X41);
      pp.query("X41", params.unit_test.X41);

      unit_test_rp::X42 = 0.0e0_rt;
      pp.query("X42", unit_test_rp::X42);
      pp.query("X42", params.unit_test.X42);

      unit_test_rp::X43 = 0.0e0_rt;
      pp.query("X43", unit_test_rp::X43);
      pp.query("X43", params.unit_test.X43);

      unit_test_rp::X44 = 0.0e0_rt;
      pp.query("X44", unit_test_rp::X44);
      pp.query("X44", params.unit_test.X44);

      unit_test_rp::X45 = 0.0e0_rt;
      pp.query("X45", unit_test_rp::X45);
      pp.query("X45", params.unit_test.X45);

      unit_test_rp::X46 = 0.0e0_rt;
      pp.query("X46", unit_test_rp::X46);
      pp.query("X46", params.unit_test.X46);

      unit_test_rp::X47 = 0.0e0_rt;
      pp.query("X47", unit_test_rp::X47);
      pp.query("X47", params.unit_test.X47);

      unit_test_rp::X48 = 0.0e0_rt;
      pp.query("X48", unit_test_rp::X48);
      pp.query("X48", params.unit_test.X48);

      unit_test_rp::X49 = 0.0e0_rt;
      pp.query("X49", unit_test_rp::X49);
      pp.query("X49", params.unit_test.X49);

      unit_test_rp::X50 = 0.0e0_rt;
      pp.query("X50", unit_test_rp::X50);
      pp.query("X50", params.unit_test.X50);

      unit_test_rp::X51 = 0.0e0_rt;
      pp.query("X51", unit_test_rp::X51);
      pp.query("X51", params.unit_test.X51);

      unit_test_rp::X52 = 0.0e0_rt;
      pp.query("X52", unit_test_rp::X52);
      pp.query("X52", params.unit_test.X52);

      unit_test_rp::X53 = 0.0e0_rt;
      pp.query("X53", unit_test_rp::X53);
      pp.query("X53", params.unit_test.X53);

      unit_test_rp::X54 = 0.0e0_rt;
      pp.query("X54", unit_test_rp::X54);
      pp.query("X54", params.unit_test.X54);

      unit_test_rp::X55 = 0.0e0_rt;
      pp.query("X55", unit_test_rp::X55);
      pp.query("X55", params.unit_test.X55);

      unit_test_rp::X56 = 0.0e0_rt;
      pp.query("X56", unit_test_rp::X56);
      pp.query("X56", params.unit_test.X56);

      unit_test_rp::X57 = 0.0e0_rt;
      pp.query("X57", unit_test_rp::X57);
      pp.query("X57", params.unit_test.X57);

      unit_test_rp::X58 = 0.0e0_rt;
      pp.query("X58", unit_test_rp::X58);
      pp.query("X58", params.unit_test.X58);

      unit_test_rp::X59 = 0.0e0_rt;
      pp.query("X59", unit_test_rp::X59);
      pp.query("X59", params.unit_test.X59);

      unit_test_rp::X60 = 0.0e0_rt;
      pp.query("X60", unit_test_rp::X60);
      pp.query("X60", params.unit_test.X60);

      unit_test_rp::X61 = 0.0e0_rt;
      pp.query("X61", unit_test_rp::X61);
      pp.query("X61", params.unit_test.X61);

      unit_test_rp::X62 = 0.0e0_rt;
      pp.query("X62", unit_test_rp::X62);
      pp.query("X62", params.unit_test.X62);

      unit_test_rp::X63 = 0.0e0_rt;
      pp.query("X63", unit_test_rp::X63);
      pp.query("X63", params.unit_test.X63);

      unit_test_rp::X64 = 0.0e0_rt;
      pp.query("X64", unit_test_rp::X64);
      pp.query("X64", params.unit_test.X64);

      unit_test_rp::X65 = 0.0e0_rt;
      pp.query("X65", unit_test_rp::X65);
      pp.query("X65", params.unit_test.X65);

      unit_test_rp::X66 = 0.0e0_rt;
      pp.query("X66", unit_test_rp::X66);
      pp.query("X66", params.unit_test.X66);

      unit_test_rp::X67 = 0.0e0_rt;
      pp.query("X67", unit_test_rp::X67);
      pp.query("X67", params.unit_test.X67);

      unit_test_rp::X68 = 0.0e0_rt;
      pp.query("X68", unit_test_rp::X68);
      pp.query("X68", params.unit_test.X68);

      unit_test_rp::X69 = 0.0e0_rt;
      pp.query("X69", unit_test_rp::X69);
      pp.query("X69", params.unit_test.X69);

      unit_test_rp::X70 = 0.0e0_rt;
      pp.query("X70", unit_test_rp::X70);
      pp.query("X70", params.unit_test.X70);

      unit_test_rp::X71 = 0.0e0_rt;
      pp.query("X71", unit_test_rp::X71);
      pp.query("X71", params.unit_test.X71);

      unit_test_rp::X72 = 0.0e0_rt;
      pp.query("X72", unit_test_rp::X72);
      pp.query("X72", params.unit_test.X72);

      unit_test_rp::X73 = 0.0e0_rt;
      pp.query("X73", unit_test_rp::X73);
      pp.query("X73", params.unit_test.X73);

      unit_test_rp::X74 = 0.0e0_rt;
      pp.query("X74", unit_test_rp::X74);
      pp.query("X74", params.unit_test.X74);

      unit_test_rp::X75 = 0.0e0_rt;
      pp.query("X75", unit_test_rp::X75);
      pp.query("X75", params.unit_test.X75);

      unit_test_rp::X76 = 0.0e0_rt;
      pp.query("X76", unit_test_rp::X76);
      pp.query("X76", params.unit_test.X76);

      unit_test_rp::X77 = 0.0e0_rt;
      pp.query("X77", unit_test_rp::X77);
      pp.query("X77", params.unit_test.X77);

      unit_test_rp::X78 = 0.0e0_rt;
      pp.query("X78", unit_test_rp::X78);
      pp.query("X78", params.unit_test.X78);

      unit_test_rp::X79 = 0.0e0_rt;
      pp.query("X79", unit_test_rp::X79);
      pp.query("X79", params.unit_test.X79);

      unit_test_rp::X80 = 0.0e0_rt;
      pp.query("X80", unit_test_rp::X80);
      pp.query("X80", params.unit_test.X80);

      unit_test_rp::X81 = 0.0e0_rt;
      pp.query("X81", unit_test_rp::X81);
      pp.query("X81", params.unit_test.X81);

      unit_test_rp::X82 = 0.0e0_rt;
      pp.query("X82", unit_test_rp::X82);
      pp.query("X82", params.unit_test.X82);

      unit_test_rp::X83 = 0.0e0_rt;
      pp.query("X83", unit_test_rp::X83);
      pp.query("X83", params.unit_test.X83);

      unit_test_rp::X84 = 0.0e0_rt;
      pp.query("X84", unit_test_rp::X84);
      pp.query("X84", params.unit_test.X84);

      unit_test_rp::X85 = 0.0e0_rt;
      pp.query("X85", unit_test_rp::X85);
      pp.query("X85", params.unit_test.X85);

      unit_test_rp::X86 = 0.0e0_rt;
      pp.query("X86", unit_test_rp::X86);
      pp.query("X86", params.unit_test.X86);

      unit_test_rp::X87 = 0.0e0_rt;
      pp.query("X87", unit_test_rp::X87);
      pp.query("X87", params.unit_test.X87);

      unit_test_rp::X88 = 0.0e0_rt;
      pp.query("X88", unit_test_rp::X88);
      pp.query("X88", params.unit_test.X88);

      unit_test_rp::X89 = 0.0e0_rt;
      pp.query("X89", unit_test_rp::X89);
      pp.query("X89", params.unit_test.X89);

      unit_test_rp::X90 = 0.0e0_rt;
      pp.query("X90", unit_test_rp::X90);
      pp.query("X90", params.unit_test.X90);

      unit_test_rp::X91 = 0.0e0_rt;
      pp.query("X91", unit_test_rp::X91);
      pp.query("X91", params.unit_test.X91);

      unit_test_rp::X92 = 0.0e0_rt;
      pp.query("X92", unit_test_rp::X92);
      pp.query("X92", params.unit_test.X92);

      unit_test_rp::X93 = 0.0e0_rt;
      pp.query("X93", unit_test_rp::X93);
      pp.query("X93", params.unit_test.X93);

      unit_test_rp::X94 = 0.0e0_rt;
      pp.query("X94", unit_test_rp::X94);
      pp.query("X94", params.unit_test.X94);

      unit_test_rp::X95 = 0.0e0_rt;
      pp.query("X95", unit_test_rp::X95);
      pp.query("X95", params.unit_test.X95);

      unit_test_rp::X96 = 0.0e0_rt;
      pp.query("X96", unit_test_rp::X96);
      pp.query("X96", params.unit_test.X96);

      unit_test_rp::X97 = 0.0e0_rt;
      pp.query("X97", unit_test_rp::X97);
      pp.query("X97", params.unit_test.X97);

      unit_test_rp::X98 = 0.0e0_rt;
      pp.query("X98", unit_test_rp::X98);
      pp.query("X98", params.unit_test.X98);

      unit_test_rp::X99 = 0.0e0_rt;
      pp.query("X99", unit_test_rp::X99);
      pp.query("X99", params.unit_test.X99);

      unit_test_rp::X100 = 0.0e0_rt;
      pp.query("X100", unit_test_rp::X100);
      pp.query("X100", params.unit_test.X100);

      unit_test_rp::X101 = 0.0e0_rt;
      pp.query("X101", unit_test_rp::X101);
      pp.query("X101", params.unit_test.X101);

      unit_test_rp::X102 = 0.0e0_rt;
      pp.query("X102", unit_test_rp::X102);
      pp.query("X102", params.unit_test.X102);

      unit_test_rp::X103 = 0.0e0_rt;
      pp.query("X103", unit_test_rp::X103);
      pp.query("X103", params.unit_test.X103);

      unit_test_rp::X104 = 0.0e0_rt;
      pp.query("X104", unit_test_rp::X104);
      pp.query("X104", params.unit_test.X104);

      unit_test_rp::X105 = 0.0e0_rt;
      pp.query("X105", unit_test_rp::X105);
      pp.query("X105", params.unit_test.X105);

      unit_test_rp::X106 = 0.0e0_rt;
      pp.query("X106", unit_test_rp::X106);
      pp.query("X106", params.unit_test.X106);

      unit_test_rp::X107 = 0.0e0_rt;
      pp.query("X107", unit_test_rp::X107);
      pp.query("X107", params.unit_test.X107);

      unit_test_rp::X108 = 0.0e0_rt;
      pp.query("X108", unit_test_rp::X108);
      pp.query("X108", params.unit_test.X108);

      unit_test_rp::X109 = 0.0e0_rt;
      pp.query("X109", unit_test_rp::X109);
      pp.query("X109", params.unit_test.X109);

      unit_test_rp::X110 = 0.0e0_rt;
      pp.query("X110", unit_test_rp::X110);
      pp.query("X110", params.unit_test.X110);

      unit_test_rp::X111 = 0.0e0_rt;
      pp.query("X111", unit_test_rp::X111);
      pp.query("X111", params.unit_test.X111);

      unit_test_rp::X112 = 0.0e0_rt;
      pp.query("X112", unit_test_rp::X112);
      pp.query("X112", params.unit_test.X112);

      unit_test_rp::X113 = 0.0e0_rt;
      pp.query("X113", unit_test_rp::X113);
      pp.query("X113", params.unit_test.X113);

      unit_test_rp::X114 = 0.0e0_rt;
      pp.query("X114", unit_test_rp::X114);
      pp.query("X114", params.unit_test.X114);

      unit_test_rp::X115 = 0.0e0_rt;
      pp.query("X115", unit_test_rp::X115);
      pp.query("X115", params.unit_test.X115);

      unit_test_rp::X116 = 0.0e0_rt;
      pp.query("X116", unit_test_rp::X116);
      pp.query("X116", params.unit_test.X116);

      unit_test_rp::X117 = 0.0e0_rt;
      pp.query("X117", unit_test_rp::X117);
      pp.query("X117", params.unit_test.X117);

      unit_test_rp::X118 = 0.0e0_rt;
      pp.query("X118", unit_test_rp::X118);
      pp.query("X118", params.unit_test.X118);

      unit_test_rp::X119 = 0.0e0_rt;
      pp.query("X119", unit_test_rp::X119);
      pp.query("X119", params.unit_test.X119);

      unit_test_rp::X120 = 0.0e0_rt;
      pp.query("X120", unit_test_rp::X120);
      pp.query("X120", params.unit_test.X120);

      unit_test_rp::X121 = 0.0e0_rt;
      pp.query("X121", unit_test_rp::X121);
      pp.query("X121", params.unit_test.X121);

      unit_test_rp::X122 = 0.0e0_rt;
      pp.query("X122", unit_test_rp::X122);
      pp.query("X122", params.unit_test.X122);

      unit_test_rp::X123 = 0.0e0_rt;
      pp.query("X123", unit_test_rp::X123);
      pp.query("X123", params.unit_test.X123);

      unit_test_rp::X124 = 0.0e0_rt;
      pp.query("X124", unit_test_rp::X124);
      pp.query("X124", params.unit_test.X124);

      unit_test_rp::X125 = 0.0e0_rt;
      pp.query("X125", unit_test_rp::X125);
      pp.query("X125", params.unit_test.X125);

      unit_test_rp::X126 = 0.0e0_rt;
      pp.query("X126", unit_test_rp::X126);
      pp.query("X126", params.unit_test.X126);

      unit_test_rp::X127 = 0.0e0_rt;
      pp.query("X127", unit_test_rp::X127);
      pp.query("X127", params.unit_test.X127);

      unit_test_rp::X128 = 0.0e0_rt;
      pp.query("X128", unit_test_rp::X128);
      pp.query("X128", params.unit_test.X128);

      unit_test_rp::X129 = 0.0e0_rt;
      pp.query("X129", unit_test_rp::X129);
      pp.query("X129", params.unit_test.X129);

      unit_test_rp::X130 = 0.0e0_rt;
      pp.query("X130", unit_test_rp::X130);
      pp.query("X130", params.unit_test.X130);

      unit_test_rp::X131 = 0.0e0_rt;
      pp.query("X131", unit_test_rp::X131);
      pp.query("X131", params.unit_test.X131);

      unit_test_rp::X132 = 0.0e0_rt;
      pp.query("X132", unit_test_rp::X132);
      pp.query("X132", params.unit_test.X132);

      unit_test_rp::X133 = 0.0e0_rt;
      pp.query("X133", unit_test_rp::X133);
      pp.query("X133", params.unit_test.X133);

      unit_test_rp::X134 = 0.0e0_rt;
      pp.query("X134", unit_test_rp::X134);
      pp.query("X134", params.unit_test.X134);

      unit_test_rp::X135 = 0.0e0_rt;
      pp.query("X135", unit_test_rp::X135);
      pp.query("X135", params.unit_test.X135);

      unit_test_rp::X136 = 0.0e0_rt;
      pp.query("X136", unit_test_rp::X136);
      pp.query("X136", params.unit_test.X136);

      unit_test_rp::X137 = 0.0e0_rt;
      pp.query("X137", unit_test_rp::X137);
      pp.query("X137", params.unit_test.X137);

      unit_test_rp::X138 = 0.0e0_rt;
      pp.query("X138", unit_test_rp::X138);
      pp.query("X138", params.unit_test.X138);

      unit_test_rp::X139 = 0.0e0_rt;
      pp.query("X139", unit_test_rp::X139);
      pp.query("X139", params.unit_test.X139);

      unit_test_rp::X140 = 0.0e0_rt;
      pp.query("X140", unit_test_rp::X140);
      pp.query("X140", params.unit_test.X140);

      unit_test_rp::X141 = 0.0e0_rt;
      pp.query("X141", unit_test_rp::X141);
      pp.query("X141", params.unit_test.X141);

      unit_test_rp::X142 = 0.0e0_rt;
      pp.query("X142", unit_test_rp::X142);
      pp.query("X142", params.unit_test.X142);

      unit_test_rp::X143 = 0.0e0_rt;
      pp.query("X143", unit_test_rp::X143);
      pp.query("X143", params.unit_test.X143);

      unit_test_rp::X144 = 0.0e0_rt;
      pp.query("X144", unit_test_rp::X144);
      pp.query("X144", params.unit_test.X144);

      unit_test_rp::X145 = 0.0e0_rt;
      pp.query("X145", unit_test_rp::X145);
      pp.query("X145", params.unit_test.X145);

      unit_test_rp::X146 = 0.0e0_rt;
      pp.query("X146", unit_test_rp::X146);
      pp.query("X146", params.unit_test.X146);

      unit_test_rp::X147 = 0.0e0_rt;
      pp.query("X147", unit_test_rp::X147);
      pp.query("X147", params.unit_test.X147);

      unit_test_rp::X148 = 0.0e0_rt;
      pp.query("X148", unit_test_rp::X148);
      pp.query("X148", params.unit_test.X148);

      unit_test_rp::X149 = 0.0e0_rt;
      pp.query("X149", unit_test_rp::X149);
      pp.query("X149", params.unit_test.X149);

      unit_test_rp::X150 = 0.0e0_rt;
      pp.query("X150", unit_test_rp::X150);
      pp.query("X150", params.unit_test.X150);

      unit_test_rp::X151 = 0.0e0_rt;
      pp.query("X151", unit_test_rp::X151);
      pp.query("X151", params.unit_test.X151);

      unit_test_rp::X152 = 0.0e0_rt;
      pp.query("X152", unit_test_rp::X152);
      pp.query("X152", params.unit_test.X152);

      unit_test_rp::X153 = 0.0e0_rt;
      pp.query("X153", unit_test_rp::X153);
      pp.query("X153", params.unit_test.X153);

      unit_test_rp::X154 = 0.0e0_rt;
      pp.query("X154", unit_test_rp::X154);
      pp.query("X154", params.unit_test.X154);

      unit_test_rp::X155 = 0.0e0_rt;
      pp.query("X155", unit_test_rp::X155);
      pp.query("X155", params.unit_test.X155);

      unit_test_rp::X156 = 0.0e0_rt;
      pp.query("X156", unit_test_rp::X156);
      pp.query("X156", params.unit_test.X156);

      unit_test_rp::X157 = 0.0e0_rt;
      pp.query("X157", unit_test_rp::X157);
      pp.query("X157", params.unit_test.X157);

      unit_test_rp::X158 = 0.0e0_rt;
      pp.query("X158", unit_test_rp::X158);
      pp.query("X158", params.unit_test.X158);

      unit_test_rp::X159 = 0.0e0_rt;
      pp.query("X159", unit_test_rp::X159);
      pp.query("X159", params.unit_test.X159);

      unit_test_rp::X160 = 0.0e0_rt;
      pp.query("X160", unit_test_rp::X160);
      pp.query("X160", params.unit_test.X160);

      unit_test_rp::X161 = 0.0e0_rt;
      pp.query("X161", unit_test_rp::X161);
      pp.query("X161", params.unit_test.X161);

      unit_test_rp::X162 = 0.0e0_rt;
      pp.query("X162", unit_test_rp::X162);
      pp.query("X162", params.unit_test.X162);

      unit_test_rp::X163 = 0.0e0_rt;
      pp.query("X163", unit_test_rp::X163);
      pp.query("X163", params.unit_test.X163);

      unit_test_rp::X164 = 0.0e0_rt;
      pp.query("X164", unit_test_rp::X164);
      pp.query("X164", params.unit_test.X164);

      unit_test_rp::X165 = 0.0e0_rt;
      pp.query("X165", unit_test_rp::X165);
      pp.query("X165", params.unit_test.X165);

      unit_test_rp::X166 = 0.0e0_rt;
      pp.query("X166", unit_test_rp::X166);
      pp.query("X166", params.unit_test.X166);

      unit_test_rp::X167 = 0.0e0_rt;
      pp.query("X167", unit_test_rp::X167);
      pp.query("X167", params.unit_test.X167);

      unit_test_rp::X168 = 0.0e0_rt;
      pp.query("X168", unit_test_rp::X168);
      pp.query("X168", params.unit_test.X168);

      unit_test_rp::X169 = 0.0e0_rt;
      pp.query("X169", unit_test_rp::X169);
      pp.query("X169", params.unit_test.X169);

      unit_test_rp::X170 = 0.0e0_rt;
      pp.query("X170", unit_test_rp::X170);
      pp.query("X170", params.unit_test.X170);

      unit_test_rp::X171 = 0.0e0_rt;
      pp.query("X171", unit_test_rp::X171);
      pp.query("X171", params.unit_test.X171);

      unit_test_rp::X172 = 0.0e0_rt;
      pp.query("X172", unit_test_rp::X172);
      pp.query("X172", params.unit_test.X172);

      unit_test_rp::X173 = 0.0e0_rt;
      pp.query("X173", unit_test_rp::X173);
      pp.query("X173", params.unit_test.X173);

      unit_test_rp::X174 = 0.0e0_rt;
      pp.query("X174", unit_test_rp::X174);
      pp.query("X174", params.unit_test.X174);

      unit_test_rp::X175 = 0.0e0_rt;
      pp.query("X175", unit_test_rp::X175);
      pp.query("X175", params.unit_test.X175);

      unit_test_rp::X176 = 0.0e0_rt;
      pp.query("X176", unit_test_rp::X176);
      pp.query("X176", params.unit_test.X176);

      unit_test_rp::X177 = 0.0e0_rt;
      pp.query("X177", unit_test_rp::X177);
      pp.query("X177", params.unit_test.X177);

      unit_test_rp::X178 = 0.0e0_rt;
      pp.query("X178", unit_test_rp::X178);
      pp.query("X178", params.unit_test.X178);

      unit_test_rp::X179 = 0.0e0_rt;
      pp.query("X179", unit_test_rp::X179);
      pp.query("X179", params.unit_test.X179);

      unit_test_rp::X180 = 0.0e0_rt;
      pp.query("X180", unit_test_rp::X180);
      pp.query("X180", params.unit_test.X180);

      unit_test_rp::X181 = 0.0e0_rt;
      pp.query("X181", unit_test_rp::X181);
      pp.query("X181", params.unit_test.X181);

      unit_test_rp::X182 = 0.0e0_rt;
      pp.query("X182", unit_test_rp::X182);
      pp.query("X182", params.unit_test.X182);

      unit_test_rp::X183 = 0.0e0_rt;
      pp.query("X183", unit_test_rp::X183);
      pp.query("X183", params.unit_test.X183);

      unit_test_rp::X184 = 0.0e0_rt;
      pp.query("X184", unit_test_rp::X184);
      pp.query("X184", params.unit_test.X184);

      unit_test_rp::X185 = 0.0e0_rt;
      pp.query("X185", unit_test_rp::X185);
      pp.query("X185", params.unit_test.X185);

      unit_test_rp::X186 = 0.0e0_rt;
      pp.query("X186", unit_test_rp::X186);
      pp.query("X186", params.unit_test.X186);

      unit_test_rp::X187 = 0.0e0_rt;
      pp.query("X187", unit_test_rp::X187);
      pp.query("X187", params.unit_test.X187);

      unit_test_rp::X188 = 0.0e0_rt;
      pp.query("X188", unit_test_rp::X188);
      pp.query("X188", params.unit_test.X188);

      unit_test_rp::X189 = 0.0e0_rt;
      pp.query("X189", unit_test_rp::X189);
      pp.query("X189", params.unit_test.X189);

      unit_test_rp::X190 = 0.0e0_rt;
      pp.query("X190", unit_test_rp::X190);
      pp.query("X190", params.unit_test.X190);

      unit_test_rp::X191 = 0.0e0_rt;
      pp.query("X191", unit_test_rp::X191);
      pp.query("X191", params.unit_test.X191);

      unit_test_rp::X192 = 0.0e0_rt;
      pp.query("X192", unit_test_rp::X192);
      pp.query("X192", params.unit_test.X192);

      unit_test_rp::X193 = 0.0e0_rt;
      pp.query("X193", unit_test_rp::X193);
      pp.query("X193", params.unit_test.X193);

      unit_test_rp::X194 = 0.0e0_rt;
      pp.query("X194", unit_test_rp::X194);
      pp.query("X194", params.unit_test.X194);

      unit_test_rp::X195 = 0.0e0_rt;
      pp.query("X195", unit_test_rp::X195);
      pp.query("X195", params.unit_test.X195);

      unit_test_rp::X196 = 0.0e0_rt;
      pp.query("X196", unit_test_rp::X196);
      pp.query("X196", params.unit_test.X196);

      unit_test_rp::X197 = 0.0e0_rt;
      pp.query("X197", unit_test_rp::X197);
      pp.query("X197", params.unit_test.X197);

      unit_test_rp::X198 = 0.0e0_rt;
      pp.query("X198", unit_test_rp::X198);
      pp.query("X198", params.unit_test.X198);

      unit_test_rp::X199 = 0.0e0_rt;
      pp.query("X199", unit_test_rp::X199);
      pp.query("X199", params.unit_test.X199);

      unit_test_rp::X200 = 0.0e0_rt;
      pp.query("X200", unit_test_rp::X200);
      pp.query("X200", params.unit_test.X200);

      unit_test_rp::X201 = 0.0e0_rt;
      pp.query("X201", unit_test_rp::X201);
      pp.query("X201", params.unit_test.X201);

      unit_test_rp::X202 = 0.0e0_rt;
      pp.query("X202", unit_test_rp::X202);
      pp.query("X202", params.unit_test.X202);

      unit_test_rp::X203 = 0.0e0_rt;
      pp.query("X203", unit_test_rp::X203);
      pp.query("X203", params.unit_test.X203);

      unit_test_rp::X204 = 0.0e0_rt;
      pp.query("X204", unit_test_rp::X204);
      pp.query("X204", params.unit_test.X204);

      unit_test_rp::X205 = 0.0e0_rt;
      pp.query("X205", unit_test_rp::X205);
      pp.query("X205", params.unit_test.X205);

      unit_test_rp::X206 = 0.0e0_rt;
      pp.query("X206", unit_test_rp::X206);
      pp.query("X206", params.unit_test.X206);

      unit_test_rp::X207 = 0.0e0_rt;
      pp.query("X207", unit_test_rp::X207);
      pp.query("X207", params.unit_test.X207);

      unit_test_rp::X208 = 0.0e0_rt;
      pp.query("X208", unit_test_rp::X208);
      pp.query("X208", params.unit_test.X208);

      unit_test_rp::X209 = 0.0e0_rt;
      pp.query("X209", unit_test_rp::X209);
      pp.query("X209", params.unit_test.X209);

      unit_test_rp::X210 = 0.0e0_rt;
      pp.query("X210", unit_test_rp::X210);
      pp.query("X210", params.unit_test.X210);

      unit_test_rp::X211 = 0.0e0_rt;
      pp.query("X211", unit_test_rp::X211);
      pp.query("X211", params.unit_test.X211);

      unit_test_rp::X212 = 0.0e0_rt;
      pp.query("X212", unit_test_rp::X212);
      pp.query("X212", params.unit_test.X212);

      unit_test_rp::X213 = 0.0e0_rt;
      pp.query("X213", unit_test_rp::X213);
      pp.query("X213", params.unit_test.X213);

      unit_test_rp::X214 = 0.0e0_rt;
      pp.query("X214", unit_test_rp::X214);
      pp.query("X214", params.unit_test.X214);

      unit_test_rp::X215 = 0.0e0_rt;
      pp.query("X215", unit_test_rp::X215);
      pp.query("X215", params.unit_test.X215);

      unit_test_rp::X216 = 0.0e0_rt;
      pp.query("X216", unit_test_rp::X216);
      pp.query("X216", params.unit_test.X216);

      unit_test_rp::X217 = 0.0e0_rt;
      pp.query("X217", unit_test_rp::X217);
      pp.query("X217", params.unit_test.X217);

      unit_test_rp::X218 = 0.0e0_rt;
      pp.query("X218", unit_test_rp::X218);
      pp.query("X218", params.unit_test.X218);

      unit_test_rp::X219 = 0.0e0_rt;
      pp.query("X219", unit_test_rp::X219);
      pp.query("X219", params.unit_test.X219);

      unit_test_rp::X220 = 0.0e0_rt;
      pp.query("X220", unit_test_rp::X220);
      pp.query("X220", params.unit_test.X220);

      unit_test_rp::X221 = 0.0e0_rt;
      pp.query("X221", unit_test_rp::X221);
      pp.query("X221", params.unit_test.X221);

      unit_test_rp::X222 = 0.0e0_rt;
      pp.query("X222", unit_test_rp::X222);
      pp.query("X222", params.unit_test.X222);

      unit_test_rp::X223 = 0.0e0_rt;
      pp.query("X223", unit_test_rp::X223);
      pp.query("X223", params.unit_test.X223);

      unit_test_rp::X224 = 0.0e0_rt;
      pp.query("X224", unit_test_rp::X224);
      pp.query("X224", params.unit_test.X224);

      unit_test_rp::X225 = 0.0e0_rt;
      pp.query("X225", unit_test_rp::X225);
      pp.query("X225", params.unit_test.X225);

      unit_test_rp::X226 = 0.0e0_rt;
      pp.query("X226", unit_test_rp::X226);
      pp.query("X226", params.unit_test.X226);

      unit_test_rp::X227 = 0.0e0_rt;
      pp.query("X227", unit_test_rp::X227);
      pp.query("X227", params.unit_test.X227);

      unit_test_rp::X228 = 0.0e0_rt;
      pp.query("X228", unit_test_rp::X228);
      pp.query("X228", params.unit_test.X228);

      unit_test_rp::X229 = 0.0e0_rt;
      pp.query("X229", unit_test_rp::X229);
      pp.query("X229", params.unit_test.X229);

      unit_test_rp::X230 = 0.0e0_rt;
      pp.query("X230", unit_test_rp::X230);
      pp.query("X230", params.unit_test.X230);

      unit_test_rp::X231 = 0.0e0_rt;
      pp.query("X231", unit_test_rp::X231);
      pp.query("X231", params.unit_test.X231);

      unit_test_rp::X232 = 0.0e0_rt;
      pp.query("X232", unit_test_rp::X232);
      pp.query("X232", params.unit_test.X232);

      unit_test_rp::X233 = 0.0e0_rt;
      pp.query("X233", unit_test_rp::X233);
      pp.query("X233", params.unit_test.X233);

      unit_test_rp::X234 = 0.0e0_rt;
      pp.query("X234", unit_test_rp::X234);
      pp.query("X234", params.unit_test.X234);

      unit_test_rp::X235 = 0.0e0_rt;
      pp.query("X235", unit_test_rp::X235);
      pp.query("X235", params.unit_test.X235);

      unit_test_rp::X236 = 0.0e0_rt;
      pp.query("X236", unit_test_rp::X236);
      pp.query("X236", params.unit_test.X236);

      unit_test_rp::X237 = 0.0e0_rt;
      pp.query("X237", unit_test_rp::X237);
      pp.query("X237", params.unit_test.X237);

      unit_test_rp::X238 = 0.0e0_rt;
      pp.query("X238", unit_test_rp::X238);
      pp.query("X238", params.unit_test.X238);

      unit_test_rp::X239 = 0.0e0_rt;
      pp.query("X239", unit_test_rp::X239);
      pp.query("X239", params.unit_test.X239);

      unit_test_rp::X240 = 0.0e0_rt;
      pp.query("X240", unit_test_rp::X240);
      pp.query("X240", params.unit_test.X240);

      unit_test_rp::X241 = 0.0e0_rt;
      pp.query("X241", unit_test_rp::X241);
      pp.query("X241", params.unit_test.X241);

      unit_test_rp::X242 = 0.0e0_rt;
      pp.query("X242", unit_test_rp::X242);
      pp.query("X242", params.unit_test.X242);

      unit_test_rp::X243 = 0.0e0_rt;
      pp.query("X243", unit_test_rp::X243);
      pp.query("X243", params.unit_test.X243);

      unit_test_rp::X244 = 0.0e0_rt;
      pp.query("X244", unit_test_rp::X244);
      pp.query("X244", params.unit_test.X244);

      unit_test_rp::X245 = 0.0e0_rt;
      pp.query("X245", unit_test_rp::X245);
      pp.query("X245", params.unit_test.X245);

      unit_test_rp::X246 = 0.0e0_rt;
      pp.query("X246", unit_test_rp::X246);
      pp.query("X246", params.unit_test.X246);

      unit_test_rp::X247 = 0.0e0_rt;
      pp.query("X247", unit_test_rp::X247);
      pp.query("X247", params.unit_test.X247);

      unit_test_rp::X248 = 0.0e0_rt;
      pp.query("X248", unit_test_rp::X248);
      pp.query("X248", params.unit_test.X248);

      unit_test_rp::X249 = 0.0e0_rt;
      pp.query("X249", unit_test_rp::X249);
      pp.query("X249", params.unit_test.X249);

      unit_test_rp::X250 = 0.0e0_rt;
      pp.query("X250", unit_test_rp::X250);
      pp.query("X250", params.unit_test.X250);

      unit_test_rp::X251 = 0.0e0_rt;
      pp.query("X251", unit_test_rp::X251);
      pp.query("X251", params.unit_test.X251);

      unit_test_rp::X252 = 0.0e0_rt;
      pp.query("X252", unit_test_rp::X252);
      pp.query("X252", params.unit_test.X252);

      unit_test_rp::X253 = 0.0e0_rt;
      pp.query("X253", unit_test_rp::X253);
      pp.query("X253", params.unit_test.X253);

      unit_test_rp::X254 = 0.0e0_rt;
      pp.query("X254", unit_test_rp::X254);
      pp.query("X254", params.unit_test.X254);

      unit_test_rp::X255 = 0.0e0_rt;
      pp.query("X255", unit_test_rp::X255);
      pp.query("X255", params.unit_test.X255);

      unit_test_rp::X256 = 0.0e0_rt;
      pp.query("X256", unit_test_rp::X256);
      pp.query("X256", params.unit_test.X256);

      unit_test_rp::X257 = 0.0e0_rt;
      pp.query("X257", unit_test_rp::X257);
      pp.query("X257", params.unit_test.X257);

      unit_test_rp::X258 = 0.0e0_rt;
      pp.query("X258", unit_test_rp::X258);
      pp.query("X258", params.unit_test.X258);

      unit_test_rp::X259 = 0.0e0_rt;
      pp.query("X259", unit_test_rp::X259);
      pp.query("X259", params.unit_test.X259);

      unit_test_rp::X260 = 0.0e0_rt;
      pp.query("X260", unit_test_rp::X260);
      pp.query("X260", params.unit_test.X260);

      unit_test_rp::X261 = 0.0e0_rt;
      pp.query("X261", unit_test_rp::X261);
      pp.query("X261", params.unit_test.X261);

      unit_test_rp::X262 = 0.0e0_rt;
      pp.query("X262", unit_test_rp::X262);
      pp.query("X262", params.unit_test.X262);

      unit_test_rp::X263 = 0.0e0_rt;
      pp.query("X263", unit_test_rp::X263);
      pp.query("X263", params.unit_test.X263);

      unit_test_rp::X264 = 0.0e0_rt;
      pp.query("X264", unit_test_rp::X264);
      pp.query("X264", params.unit_test.X264);

      unit_test_rp::X265 = 0.0e0_rt;
      pp.query("X265", unit_test_rp::X265);
      pp.query("X265", params.unit_test.X265);

      unit_test_rp::X266 = 0.0e0_rt;
      pp.query("X266", unit_test_rp::X266);
      pp.query("X266", params.unit_test.X266);

      unit_test_rp::X267 = 0.0e0_rt;
      pp.query("X267", unit_test_rp::X267);
      pp.query("X267", params.unit_test.X267);

      unit_test_rp::X268 = 0.0e0_rt;
      pp.query("X268", unit_test_rp::X268);
      pp.query("X268", params.unit_test.X268);

      unit_test_rp::X269 = 0.0e0_rt;
      pp.query("X269", unit_test_rp::X269);
      pp.query("X269", params.unit_test.X269);

      unit_test_rp::X270 = 0.0e0_rt;
      pp.query("X270", unit_test_rp::X270);
      pp.query("X270", params.unit_test.X270);

      unit_test_rp::X271 = 0.0e0_rt;
      pp.query("X271", unit_test_rp::X271);
      pp.query("X271", params.unit_test.X271);

      unit_test_rp::X272 = 0.0e0_rt;
      pp.query("X272", unit_test_rp::X272);
      pp.query("X272", params.unit_test.X272);

      unit_test_rp::X273 = 0.0e0_rt;
      pp.query("X273", unit_test_rp::X273);
      pp.query("X273", params.unit_test.X273);

      unit_test_rp::X274 = 0.0e0_rt;
      pp.query("X274", unit_test_rp::X274);
      pp.query("X274", params.unit_test.X274);

      unit_test_rp::X275 = 0.0e0_rt;
      pp.query("X275", unit_test_rp::X275);
      pp.query("X275", params.unit_test.X275);

      unit_test_rp::X276 = 0.0e0_rt;
      pp.query("X276", unit_test_rp::X276);
      pp.query("X276", params.unit_test.X276);

      unit_test_rp::X277 = 0.0e0_rt;
      pp.query("X277", unit_test_rp::X277);
      pp.query("X277", params.unit_test.X277);

      unit_test_rp::X278 = 0.0e0_rt;
      pp.query("X278", unit_test_rp::X278);
      pp.query("X278", params.unit_test.X278);

      unit_test_rp::X279 = 0.0e0_rt;
      pp.query("X279", unit_test_rp::X279);
      pp.query("X279", params.unit_test.X279);

      unit_test_rp::X280 = 0.0e0_rt;
      pp.query("X280", unit_test_rp::X280);
      pp.query("X280", params.unit_test.X280);

      unit_test_rp::X281 = 0.0e0_rt;
      pp.query("X281", unit_test_rp::X281);
      pp.query("X281", params.unit_test.X281);

      unit_test_rp::X282 = 0.0e0_rt;
      pp.query("X282", unit_test_rp::X282);
      pp.query("X282", params.unit_test.X282);

      unit_test_rp::X283 = 0.0e0_rt;
      pp.query("X283", unit_test_rp::X283);
      pp.query("X283", params.unit_test.X283);

      unit_test_rp::X284 = 0.0e0_rt;
      pp.query("X284", unit_test_rp::X284);
      pp.query("X284", params.unit_test.X284);

      unit_test_rp::X285 = 0.0e0_rt;
      pp.query("X285", unit_test_rp::X285);
      pp.query("X285", params.unit_test.X285);

      unit_test_rp::X286 = 0.0e0_rt;
      pp.query("X286", unit_test_rp::X286);
      pp.query("X286", params.unit_test.X286);

      unit_test_rp::X287 = 0.0e0_rt;
      pp.query("X287", unit_test_rp::X287);
      pp.query("X287", params.unit_test.X287);

      unit_test_rp::X288 = 0.0e0_rt;
      pp.query("X288", unit_test_rp::X288);
      pp.query("X288", params.unit_test.X288);

      unit_test_rp::X289 = 0.0e0_rt;
      pp.query("X289", unit_test_rp::X289);
      pp.query("X289", params.unit_test.X289);

      unit_test_rp::X290 = 0.0e0_rt;
      pp.query("X290", unit_test_rp::X290);
      pp.query("X290", params.unit_test.X290);

      unit_test_rp::X291 = 0.0e0_rt;
      pp.query("X291", unit_test_rp::X291);
      pp.query("X291", params.unit_test.X291);

      unit_test_rp::X292 = 0.0e0_rt;
      pp.query("X292", unit_test_rp::X292);
      pp.query("X292", params.unit_test.X292);

      unit_test_rp::X293 = 0.0e0_rt;
      pp.query("X293", unit_test_rp::X293);
      pp.query("X293", params.unit_test.X293);

      unit_test_rp::X294 = 0.0e0_rt;
      pp.query("X294", unit_test_rp::X294);
      pp.query("X294", params.unit_test.X294);

      unit_test_rp::X295 = 0.0e0_rt;
      pp.query("X295", unit_test_rp::X295);
      pp.query("X295", params.unit_test.X295);

      unit_test_rp::X296 = 0.0e0_rt;
      pp.query("X296", unit_test_rp::X296);
      pp.query("X296", params.unit_test.X296);

      unit_test_rp::X297 = 0.0e0_rt;
      pp.query("X297", unit_test_rp::X297);
      pp.query("X297", params.unit_test.X297);

      unit_test_rp::X298 = 0.0e0_rt;
      pp.query("X298", unit_test_rp::X298);
      pp.query("X298", params.unit_test.X298);

      unit_test_rp::X299 = 0.0e0_rt;
      pp.query("X299", unit_test_rp::X299);
      pp.query("X299", params.unit_test.X299);

      unit_test_rp::X300 = 0.0e0_rt;
      pp.query("X300", unit_test_rp::X300);
      pp.query("X300", params.unit_test.X300);

      unit_test_rp::X301 = 0.0e0_rt;
      pp.query("X301", unit_test_rp::X301);
      pp.query("X301", params.unit_test.X301);

      unit_test_rp::X302 = 0.0e0_rt;
      pp.query("X302", unit_test_rp::X302);
      pp.query("X302", params.unit_test.X302);

      unit_test_rp::X303 = 0.0e0_rt;
      pp.query("X303", unit_test_rp::X303);
      pp.query("X303", params.unit_test.X303);

      unit_test_rp::X304 = 0.0e0_rt;
      pp.query("X304", unit_test_rp::X304);
      pp.query("X304", params.unit_test.X304);

      unit_test_rp::X305 = 0.0e0_rt;
      pp.query("X305", unit_test_rp::X305);
      pp.query("X305", params.unit_test.X305);

      unit_test_rp::X306 = 0.0e0_rt;
      pp.query("X306", unit_test_rp::X306);
      pp.query("X306", params.unit_test.X306);

      unit_test_rp::X307 = 0.0e0_rt;
      pp.query("X307", unit_test_rp::X307);
      pp.query("X307", params.unit_test.X307);

      unit_test_rp::X308 = 0.0e0_rt;
      pp.query("X308", unit_test_rp::X308);
      pp.query("X308", params.unit_test.X308);

      unit_test_rp::X309 = 0.0e0_rt;
      pp.query("X309", unit_test_rp::X309);
      pp.query("X309", params.unit_test.X309);

      unit_test_rp::X310 = 0.0e0_rt;
      pp.query("X310", unit_test_rp::X310);
      pp.query("X310", params.unit_test.X310);

      unit_test_rp::X311 = 0.0e0_rt;
      pp.query("X311", unit_test_rp::X311);
      pp.query("X311", params.unit_test.X311);

      unit_test_rp::X312 = 0.0e0_rt;
      pp.query("X312", unit_test_rp::X312);
      pp.query("X312", params.unit_test.X312);

      unit_test_rp::X313 = 0.0e0_rt;
      pp.query("X313", unit_test_rp::X313);
      pp.query("X313", params.unit_test.X313);

      unit_test_rp::X314 = 0.0e0_rt;
      pp.query("X314", unit_test_rp::X314);
      pp.query("X314", params.unit_test.X314);

      unit_test_rp::X315 = 0.0e0_rt;
      pp.query("X315", unit_test_rp::X315);
      pp.query("X315", params.unit_test.X315);

      unit_test_rp::X316 = 0.0e0_rt;
      pp.query("X316", unit_test_rp::X316);
      pp.query("X316", params.unit_test.X316);

      unit_test_rp::X317 = 0.0e0_rt;
      pp.query("X317", unit_test_rp::X317);
      pp.query("X317", params.unit_test.X317);

      unit_test_rp::X318 = 0.0e0_rt;
      pp.query("X318", unit_test_rp::X318);
      pp.query("X318", params.unit_test.X318);

      unit_test_rp::X319 = 0.0e0_rt;
      pp.query("X319", unit_test_rp::X319);
      pp.query("X319", params.unit_test.X319);

      unit_test_rp::X320 = 0.0e0_rt;
      pp.query("X320", unit_test_rp::X320);
      pp.query("X320", params.unit_test.X320);

      unit_test_rp::X321 = 0.0e0_rt;
      pp.query("X321", unit_test_rp::X321);
      pp.query("X321", params.unit_test.X321);

      unit_test_rp::X322 = 0.0e0_rt;
      pp.query("X322", unit_test_rp::X322);
      pp.query("X322", params.unit_test.X322);

      unit_test_rp::X323 = 0.0e0_rt;
      pp.query("X323", unit_test_rp::X323);
      pp.query("X323", params.unit_test.X323);

      unit_test_rp::X324 = 0.0e0_rt;
      pp.query("X324", unit_test_rp::X324);
      pp.query("X324", params.unit_test.X324);

      unit_test_rp::X325 = 0.0e0_rt;
      pp.query("X325", unit_test_rp::X325);
      pp.query("X325", params.unit_test.X325);

      unit_test_rp::X326 = 0.0e0_rt;
      pp.query("X326", unit_test_rp::X326);
      pp.query("X326", params.unit_test.X326);

      unit_test_rp::X327 = 0.0e0_rt;
      pp.query("X327", unit_test_rp::X327);
      pp.query("X327", params.unit_test.X327);

      unit_test_rp::X328 = 0.0e0_rt;
      pp.query("X328", unit_test_rp::X328);
      pp.query("X328", params.unit_test.X328);

      unit_test_rp::X329 = 0.0e0_rt;
      pp.query("X329", unit_test_rp::X329);
      pp.query("X329", params.unit_test.X329);

      unit_test_rp::X330 = 0.0e0_rt;
      pp.query("X330", unit_test_rp::X330);
      pp.query("X330", params.unit_test.X330);

      unit_test_rp::X331 = 0.0e0_rt;
      pp.query("X331", unit_test_rp::X331);
      pp.query("X331", params.unit_test.X331);

      unit_test_rp::X332 = 0.0e0_rt;
      pp.query("X332", unit_test_rp::X332);
      pp.query("X332", params.unit_test.X332);

      unit_test_rp::X333 = 0.0e0_rt;
      pp.query("X333", unit_test_rp::X333);
      pp.query("X333", params.unit_test.X333);

      unit_test_rp::X334 = 0.0e0_rt;
      pp.query("X334", unit_test_rp::X334);
      pp.query("X334", params.unit_test.X334);

      unit_test_rp::X335 = 0.0e0_rt;
      pp.query("X335", unit_test_rp::X335);
      pp.query("X335", params.unit_test.X335);

      unit_test_rp::X336 = 0.0e0_rt;
      pp.query("X336", unit_test_rp::X336);
      pp.query("X336", params.unit_test.X336);

      unit_test_rp::X337 = 0.0e0_rt;
      pp.query("X337", unit_test_rp::X337);
      pp.query("X337", params.unit_test.X337);

      unit_test_rp::X338 = 0.0e0_rt;
      pp.query("X338", unit_test_rp::X338);
      pp.query("X338", params.unit_test.X338);

      unit_test_rp::X339 = 0.0e0_rt;
      pp.query("X339", unit_test_rp::X339);
      pp.query("X339", params.unit_test.X339);

      unit_test_rp::X340 = 0.0e0_rt;
      pp.query("X340", unit_test_rp::X340);
      pp.query("X340", params.unit_test.X340);

      unit_test_rp::X341 = 0.0e0_rt;
      pp.query("X341", unit_test_rp::X341);
      pp.query("X341", params.unit_test.X341);

      unit_test_rp::X342 = 0.0e0_rt;
      pp.query("X342", unit_test_rp::X342);
      pp.query("X342", params.unit_test.X342);

      unit_test_rp::X343 = 0.0e0_rt;
      pp.query("X343", unit_test_rp::X343);
      pp.query("X343", params.unit_test.X343);

      unit_test_rp::X344 = 0.0e0_rt;
      pp.query("X344", unit_test_rp::X344);
      pp.query("X344", params.unit_test.X344);

      unit_test_rp::X345 = 0.0e0_rt;
      pp.query("X345", unit_test_rp::X345);
      pp.query("X345", params.unit_test.X345);

      unit_test_rp::X346 = 0.0e0_rt;
      pp.query("X346", unit_test_rp::X346);
      pp.query("X346", params.unit_test.X346);

      unit_test_rp::X347 = 0.0e0_rt;
      pp.query("X347", unit_test_rp::X347);
      pp.query("X347", params.unit_test.X347);

      unit_test_rp::X348 = 0.0e0_rt;
      pp.query("X348", unit_test_rp::X348);
      pp.query("X348", params.unit_test.X348);

      unit_test_rp::X349 = 0.0e0_rt;
      pp.query("X349", unit_test_rp::X349);
      pp.query("X349", params.unit_test.X349);

      unit_test_rp::X350 = 0.0e0_rt;
      pp.query("X350", unit_test_rp::X350);
      pp.query("X350", params.unit_test.X350);

      unit_test_rp::X351 = 0.0e0_rt;
      pp.query("X351", unit_test_rp::X351);
      pp.query("X351", params.unit_test.X351);

      unit_test_rp::X352 = 0.0e0_rt;
      pp.query("X352", unit_test_rp::X352);
      pp.query("X352", params.unit_test.X352);

      unit_test_rp::X353 = 0.0e0_rt;
      pp.query("X353", unit_test_rp::X353);
      pp.query("X353", params.unit_test.X353);

      unit_test_rp::X354 = 0.0e0_rt;
      pp.query("X354", unit_test_rp::X354);
      pp.query("X354", params.unit_test.X354);

      unit_test_rp::X355 = 0.0e0_rt;
      pp.query("X355", unit_test_rp::X355);
      pp.query("X355", params.unit_test.X355);

      unit_test_rp::X356 = 0.0e0_rt;
      pp.query("X356", unit_test_rp::X356);
      pp.query("X356", params.unit_test.X356);

      unit_test_rp::X357 = 0.0e0_rt;
      pp.query("X357", unit_test_rp::X357);
      pp.query("X357", params.unit_test.X357);

      unit_test_rp::X358 = 0.0e0_rt;
      pp.query("X358", unit_test_rp::X358);
      pp.query("X358", params.unit_test.X358);

      unit_test_rp::X359 = 0.0e0_rt;
      pp.query("X359", unit_test_rp::X359);
      pp.query("X359", params.unit_test.X359);

      unit_test_rp::X360 = 0.0e0_rt;
      pp.query("X360", unit_test_rp::X360);
      pp.query("X360", params.unit_test.X360);

      unit_test_rp::X361 = 0.0e0_rt;
      pp.query("X361", unit_test_rp::X361);
      pp.query("X361", params.unit_test.X361);

      unit_test_rp::X362 = 0.0e0_rt;
      pp.query("X362", unit_test_rp::X362);
      pp.query("X362", params.unit_test.X362);

      unit_test_rp::X363 = 0.0e0_rt;
      pp.query("X363", unit_test_rp::X363);
      pp.query("X363", params.unit_test.X363);

      unit_test_rp::X364 = 0.0e0_rt;
      pp.query("X364", unit_test_rp::X364);
      pp.query("X364", params.unit_test.X364);

      unit_test_rp::X365 = 0.0e0_rt;
      pp.query("X365", unit_test_rp::X365);
      pp.query("X365", params.unit_test.X365);

      unit_test_rp::X366 = 0.0e0_rt;
      pp.query("X366", unit_test_rp::X366);
      pp.query("X366", params.unit_test.X366);

      unit_test_rp::X367 = 0.0e0_rt;
      pp.query("X367", unit_test_rp::X367);
      pp.query("X367", params.unit_test.X367);

      unit_test_rp::X368 = 0.0e0_rt;
      pp.query("X368", unit_test_rp::X368);
      pp.query("X368", params.unit_test.X368);

      unit_test_rp::X369 = 0.0e0_rt;
      pp.query("X369", unit_test_rp::X369);
      pp.query("X369", params.unit_test.X369);

      unit_test_rp::X370 = 0.0e0_rt;
      pp.query("X370", unit_test_rp::X370);
      pp.query("X370", params.unit_test.X370);

      unit_test_rp::X371 = 0.0e0_rt;
      pp.query("X371", unit_test_rp::X371);
      pp.query("X371", params.unit_test.X371);

      unit_test_rp::X372 = 0.0e0_rt;
      pp.query("X372", unit_test_rp::X372);
      pp.query("X372", params.unit_test.X372);

      unit_test_rp::X373 = 0.0e0_rt;
      pp.query("X373", unit_test_rp::X373);
      pp.query("X373", params.unit_test.X373);

      unit_test_rp::X374 = 0.0e0_rt;
      pp.query("X374", unit_test_rp::X374);
      pp.query("X374", params.unit_test.X374);

      unit_test_rp::X375 = 0.0e0_rt;
      pp.query("X375", unit_test_rp::X375);
      pp.query("X375", params.unit_test.X375);

      unit_test_rp::X376 = 0.0e0_rt;
      pp.query("X376", unit_test_rp::X376);
      pp.query("X376", params.unit_test.X376);

      unit_test_rp::X377 = 0.0e0_rt;
      pp.query("X377", unit_test_rp::X377);
      pp.query("X377", params.unit_test.X377);

      unit_test_rp::X378 = 0.0e0_rt;
      pp.query("X378", unit_test_rp::X378);
      pp.query("X378", params.unit_test.X378);

      unit_test_rp::X379 = 0.0e0_rt;
      pp.query("X379", unit_test_rp::X379);
      pp.query("X379", params.unit_test.X379);

      unit_test_rp::X380 = 0.0e0_rt;
      pp.query("X380", unit_test_rp::X380);
      pp.query("X380", params.unit_test.X380);

      unit_test_rp::X381 = 0.0e0_rt;
      pp.query("X381", unit_test_rp::X381);
      pp.query("X381", params.unit_test.X381);

      unit_test_rp::X382 = 0.0e0_rt;
      pp.query("X382", unit_test_rp::X382);
      pp.query("X382", params.unit_test.X382);

      unit_test_rp::X383 = 0.0e0_rt;
      pp.query("X383", unit_test_rp::X383);
      pp.query("X383", params.unit_test.X383);

      unit_test_rp::X384 = 0.0e0_rt;
      pp.query("X384", unit_test_rp::X384);
      pp.query("X384", params.unit_test.X384);

      unit_test_rp::X385 = 0.0e0_rt;
      pp.query("X385", unit_test_rp::X385);
      pp.query("X385", params.unit_test.X385);

      unit_test_rp::X386 = 0.0e0_rt;
      pp.query("X386", unit_test_rp::X386);
      pp.query("X386", params.unit_test.X386);

      unit_test_rp::X387 = 0.0e0_rt;
      pp.query("X387", unit_test_rp::X387);
      pp.query("X387", params.unit_test.X387);

      unit_test_rp::X388 = 0.0e0_rt;
      pp.query("X388", unit_test_rp::X388);
      pp.query("X388", params.unit_test.X388);

      unit_test_rp::X389 = 0.0e0_rt;
      pp.query("X389", unit_test_rp::X389);
      pp.query("X389", params.unit_test.X389);

      unit_test_rp::X390 = 0.0e0_rt;
      pp.query("X390", unit_test_rp::X390);
      pp.query("X390", params.unit_test.X390);

      unit_test_rp::X391 = 0.0e0_rt;
      pp.query("X391", unit_test_rp::X391);
      pp.query("X391", params.unit_test.X391);

      unit_test_rp::X392 = 0.0e0_rt;
      pp.query("X392", unit_test_rp::X392);
      pp.query("X392", params.unit_test.X392);

      unit_test_rp::X393 = 0.0e0_rt;
      pp.query("X393", unit_test_rp::X393);
      pp.query("X393", params.unit_test.X393);

      unit_test_rp::X394 = 0.0e0_rt;
      pp.query("X394", unit_test_rp::X394);
      pp.query("X394", params.unit_test.X394);

      unit_test_rp::X395 = 0.0e0_rt;
      pp.query("X395", unit_test_rp::X395);
      pp.query("X395", params.unit_test.X395);

      unit_test_rp::X396 = 0.0e0_rt;
      pp.query("X396", unit_test_rp::X396);
      pp.query("X396", params.unit_test.X396);

      unit_test_rp::X397 = 0.0e0_rt;
      pp.query("X397", unit_test_rp::X397);
      pp.query("X397", params.unit_test.X397);

      unit_test_rp::X398 = 0.0e0_rt;
      pp.query("X398", unit_test_rp::X398);
      pp.query("X398", params.unit_test.X398);

      unit_test_rp::X399 = 0.0e0_rt;
      pp.query("X399", unit_test_rp::X399);
      pp.query("X399", params.unit_test.X399);

      unit_test_rp::X400 = 0.0e0_rt;
      pp.query("X400", unit_test_rp::X400);
      pp.query("X400", params.unit_test.X400);

      unit_test_rp::uniform_xn = false;
      pp.query("uniform_xn", unit_test_rp::uniform_xn);
      pp.query("uniform_xn", params.unit_test.uniform_xn);

      unit_test_rp::small_temp = 1.e5_rt;
      pp.query("small_temp", unit_test_rp::small_temp);
      pp.query("small_temp", params.unit_test.small_temp);

      unit_test_rp::small_dens = 1.e5_rt;
      pp.query("small_dens", unit_test_rp::small_dens);
      pp.query("small_dens", params.unit_test.small_dens);

    }
    return params;

  }
