function der_st_dz2 = der_stima_dz2(lambda,mu,V)
%DER_STIMA_DZ2
%    DER_ST_DZ2 = DER_STIMA_DZ2(LAMBDA,MU,X1,X2,X3,X4,Y1,Y2,Y3,Y4,Z1,Z2,Z3,Z4)

%    This function was generated by the Symbolic Math Toolbox version 5.8.
%    02-Jun-2015 16:40:56
x1=V(1,1);
y1=V(1,2);
z1=V(1,3);

x2=V(2,1);
y2=V(2,2);
z2=V(2,3);

x3=V(3,1);
y3=V(3,2);
z3=V(3,3);

x4=V(4,1);
y4=V(4,2);
z4=V(4,3);

t13 = x3.*y4;
t14 = x4.*y3;
t53 = x2.*y3;
t54 = x3.*y2;
t55 = x2.*y4;
t56 = x4.*y2;
t2 = t13-t14+t53-t54-t55+t56;
t3 = x2.*z3;
t4 = x4.*z2;
t5 = x3.*z4;
t29 = x3.*z2;
t30 = x2.*z4;
t31 = x4.*z3;
t6 = t3+t4+t5-t29-t30-t31;
t7 = mu.*2.0;
t8 = lambda+t7;
t9 = y2.*z3;
t10 = y4.*z2;
t11 = y3.*z4;
t33 = y3.*z2;
t34 = y2.*z4;
t35 = y4.*z3;
t12 = t9+t10+t11-t33-t34-t35;
t15 = x1.*y2.*z3;
t16 = x2.*y3.*z1;
t17 = x3.*y1.*z2;
t18 = x1.*y4.*z2;
t19 = x2.*y1.*z4;
t20 = x4.*y2.*z1;
t21 = x1.*y3.*z4;
t22 = x3.*y4.*z1;
t23 = x4.*y1.*z3;
t24 = x2.*y4.*z3;
t25 = x3.*y2.*z4;
t26 = x4.*y3.*z2;
t36 = x1.*y3.*z2;
t37 = x2.*y1.*z3;
t38 = x3.*y2.*z1;
t39 = x1.*y2.*z4;
t40 = x2.*y4.*z1;
t41 = x4.*y1.*z2;
t42 = x1.*y4.*z3;
t43 = x3.*y1.*z4;
t44 = x4.*y3.*z1;
t45 = x2.*y3.*z4;
t46 = x3.*y4.*z2;
t47 = x4.*y2.*z3;
t27 = t15+t16+t17+t18+t19+t20+t21+t22+t23+t24+t25+t26-t36-t37-t38-t39-t40-t41-t42-t43-t44-t45-t46-t47;
t28 = y3-y4;
t32 = x3-x4;
t48 = 1.0./t27;
t49 = x1.*y3;
t50 = x4.*y1;
t57 = x3.*y1;
t58 = x1.*y4;
t51 = t13-t14+t49+t50-t57-t58;
t52 = 1.0./t27.^2;
t59 = x1.*z3;
t60 = x4.*z1;
t65 = x3.*z1;
t66 = x1.*z4;
t61 = t5-t31+t59+t60-t65-t66;
t62 = y1.*z3;
t63 = y4.*z1;
t67 = y3.*z1;
t68 = y1.*z4;
t64 = t11-t35+t62+t63-t67-t68;
t69 = x2.*z1;
t73 = x1.*z2;
t70 = t4-t30-t60+t66+t69-t73;
t71 = y2.*z1;
t76 = y1.*z2;
t72 = t10-t34-t63+t68+t71-t76;
t74 = y1-y4;
t75 = x1-x4;
t77 = x1.*y2;
t79 = x2.*y1;
t78 = t50+t55-t56-t58+t77-t79;
t80 = t3-t29-t59+t65-t69+t73;
t81 = t9-t33-t62+t67-t71+t76;
t82 = y1-y3;
t83 = x1-x3;
t84 = t49-t53+t54-t57-t77+t79;
t85 = lambda.*t6.*t28;
t86 = mu.*t6.*t28;
t87 = lambda.*t12.*t32;
t88 = mu.*t12.*t32;
t89 = t85+t86+t87+t88;
t90 = lambda.*t6.*t12;
t91 = mu.*t6.*t12;
t92 = t90+t91;
t93 = t51.*t52.*t92.*(1.0./6.0);
t94 = t93-t48.*t89.*(1.0./6.0);
t95 = t2.^2;
t96 = mu.*t95;
t97 = t12.^2;
t98 = t6.^2;
t99 = mu.*t2.*t51;
t100 = lambda.*t2.*t28;
t101 = mu.*t2.*t28;
t102 = t100+t101;
t103 = t48.*t102.*(1.0./6.0);
t104 = lambda.*t2.*t12;
t105 = mu.*t2.*t12;
t106 = t104+t105;
t107 = t103-t51.*t52.*t106.*(1.0./6.0);
t108 = lambda.*t2.*t32;
t109 = mu.*t2.*t32;
t110 = t108+t109;
t111 = lambda.*t2.*t6;
t112 = mu.*t2.*t6;
t113 = t111+t112;
t114 = t51.*t52.*t113.*(1.0./6.0);
t115 = t114-t48.*t110.*(1.0./6.0);
t116 = mu.*t6.*t32.*2.0;
t117 = mu.*t12.*t28.*2.0;
t118 = mu.*t98;
t119 = mu.*t97;
t120 = mu.*t32.*t61;
t121 = mu.*t28.*t64;
t122 = mu.*t6.*t61;
t123 = mu.*t12.*t64;
t124 = mu.*t32.*t70;
t125 = mu.*t6.*t75;
t126 = mu.*t28.*t72;
t127 = mu.*t12.*t74;
t128 = mu.*t6.*t70;
t129 = mu.*t12.*t72;
t130 = mu.*t6.*t83;
t131 = mu.*t12.*t82;
t132 = mu.*t6.*t80;
t133 = mu.*t12.*t81;
t134 = t8.*t28.*t64;
t135 = t120+t134;
t136 = t8.*t12.*t64;
t137 = t99+t122+t136;
t138 = t51.*t52.*t137.*(1.0./6.0);
t139 = t138-t48.*t135.*(1.0./6.0);
t140 = mu.*t28.*t61;
t141 = lambda.*t32.*t64;
t142 = t140+t141;
t143 = t48.*t142.*(1.0./6.0);
t144 = lambda.*t6.*t64;
t145 = mu.*t12.*t61;
t146 = t144+t145;
t147 = t143-t51.*t52.*t146.*(1.0./6.0);
t148 = lambda.*t2.*t64;
t149 = mu.*t12.*t51;
t150 = t148+t149;
t151 = t51.*t52.*t150.*(1.0./6.0);
t152 = t151-mu.*t28.*t48.*t51.*(1.0./6.0);
t153 = lambda.*t28.*t61;
t154 = mu.*t32.*t64;
t155 = t153+t154;
t156 = t48.*t155.*(1.0./6.0);
t157 = lambda.*t12.*t61;
t158 = mu.*t6.*t64;
t159 = t157+t158;
t160 = t156-t51.*t52.*t159.*(1.0./6.0);
t161 = t8.*t32.*t61;
t162 = t121+t161;
t163 = t6.*t8.*t61;
t164 = t99+t123+t163;
t165 = t51.*t52.*t164.*(1.0./6.0);
t166 = t165-t48.*t162.*(1.0./6.0);
t167 = lambda.*t2.*t61;
t168 = mu.*t6.*t51;
t169 = t167+t168;
t170 = mu.*t32.*t48.*t51.*(1.0./6.0);
t171 = t170-t51.*t52.*t169.*(1.0./6.0);
t172 = lambda.*t61.*t64;
t173 = mu.*t61.*t64;
t174 = t172+t173;
t175 = t51.*t52.*t174.*(1.0./6.0);
t176 = t51.^2;
t177 = mu.*t176;
t178 = t64.^2;
t179 = t61.^2;
t180 = lambda.*t12.*t51;
t181 = mu.*t2.*t64;
t182 = t180+t181;
t183 = t51.*t52.*t182.*(1.0./6.0);
t184 = t183-lambda.*t28.*t48.*t51.*(1.0./6.0);
t185 = lambda.*t6.*t51;
t186 = mu.*t2.*t61;
t187 = t185+t186;
t188 = lambda.*t32.*t48.*t51.*(1.0./6.0);
t189 = t188-t51.*t52.*t187.*(1.0./6.0);
t190 = t120+t121;
t191 = t2.*t8.*t51;
t192 = t122+t123+t191;
t193 = t51.*t52.*t192.*(1.0./6.0);
t194 = t193-t48.*t190.*(1.0./6.0);
t195 = lambda.*t51.*t64;
t196 = mu.*t51.*t64;
t197 = t195+t196;
t198 = lambda.*t51.*t61;
t199 = mu.*t51.*t61;
t200 = t198+t199;
t201 = t51.*t52.*t200.*(1.0./6.0);
t202 = mu.*t179;
t203 = mu.*t178;
t204 = mu.*t61.*t75;
t205 = mu.*t64.*t74;
t206 = mu.*t61.*t70;
t207 = mu.*t64.*t72;
t208 = mu.*t61.*t83;
t209 = mu.*t64.*t82;
t210 = mu.*t61.*t80;
t211 = mu.*t64.*t81;
t212 = t8.*t28.*t72;
t213 = t8.*t12.*t74;
t214 = t124+t125+t212+t213;
t215 = t8.*t12.*t72;
t264 = mu.*t2.*t78;
t216 = t51.*t52.*(t128+t215-t264).*(1.0./6.0);
t217 = t216-t48.*t214.*(1.0./6.0);
t218 = lambda.*t6.*t74;
t219 = mu.*t28.*(t4-t30-t60+t66+t69-t73);
t220 = lambda.*t32.*(t10-t34-t63+t68+t71-t76);
t221 = mu.*t12.*t75;
t222 = t218+t219+t220+t221;
t223 = t48.*t222.*(1.0./6.0);
t224 = lambda.*t6.*t72;
t225 = mu.*t12.*t70;
t226 = t224+t225;
t227 = t223-t51.*t52.*t226.*(1.0./6.0);
t228 = lambda.*t2.*t74;
t229 = t228-mu.*t28.*t78;
t230 = lambda.*t2.*t72;
t231 = t51.*t52.*(t230-mu.*t12.*t78).*(1.0./6.0);
t232 = t231-t48.*t229.*(1.0./6.0);
t233 = t8.*t64.*t74;
t234 = t204+t233;
t235 = t48.*t234.*(1.0./6.0);
t236 = t8.*t64.*t72;
t282 = mu.*t51.*t78;
t237 = t206+t236-t282;
t238 = t235-t51.*t52.*t237.*(1.0./6.0);
t239 = lambda.*t61.*t74;
t240 = mu.*t64.*t75;
t241 = t239+t240;
t242 = lambda.*t61.*t72;
t243 = mu.*t64.*t70;
t244 = t51.*t52.*(t242+t243).*(1.0./6.0);
t245 = t244-t48.*t241.*(1.0./6.0);
t246 = lambda.*t51.*t72;
t247 = t246-mu.*t64.*t78;
t248 = lambda.*t48.*t51.*t74.*(1.0./6.0);
t249 = t248-t51.*t52.*t247.*(1.0./6.0);
t250 = t4-t30-t60+t66+t69-t73;
t251 = t10-t34-t63+t68+t71-t76;
t252 = lambda.*t28.*t70;
t253 = mu.*t6.*t74;
t254 = lambda.*t12.*t75;
t255 = mu.*t32.*t72;
t256 = t48.*(t252+t253+t254+t255).*(1.0./6.0);
t257 = lambda.*t12.*t70;
t258 = mu.*t6.*t72;
t259 = t257+t258;
t260 = t256-t51.*t52.*t259.*(1.0./6.0);
t261 = t8.*t32.*t70;
t262 = t6.*t8.*t75;
t263 = t126+t127+t261+t262;
t265 = t6.*t8.*t70;
t266 = lambda.*t2.*t75;
t267 = t266-mu.*t32.*t78;
t268 = t48.*t267.*(1.0./6.0);
t269 = lambda.*t2.*t70;
t270 = t269-mu.*t6.*t78;
t271 = t268-t51.*t52.*t270.*(1.0./6.0);
t272 = mu.*t61.*t74;
t273 = lambda.*t64.*t75;
t274 = t272+t273;
t275 = lambda.*t64.*t70;
t276 = mu.*t61.*t72;
t277 = t51.*t52.*(t275+t276).*(1.0./6.0);
t278 = t277-t48.*t274.*(1.0./6.0);
t279 = t8.*t61.*t75;
t280 = t205+t279;
t281 = t48.*t280.*(1.0./6.0);
t283 = t8.*t61.*t70;
t284 = lambda.*t51.*t70;
t285 = t51.*t52.*(t284-mu.*t61.*t78).*(1.0./6.0);
t286 = t285-lambda.*t48.*t51.*t75.*(1.0./6.0);
t287 = lambda.*t70.*t74;
t288 = mu.*t70.*t74;
t289 = lambda.*t72.*t75;
t290 = mu.*t72.*t75;
t291 = t287+t288+t289+t290;
t292 = lambda.*(t4-t30-t60+t66+t69-t73).*(t10-t34-t63+t68+t71-t76);
t293 = mu.*(t4-t30-t60+t66+t69-t73).*(t10-t34-t63+t68+t71-t76);
t294 = t292+t293;
t295 = t51.*t52.*t294.*(1.0./6.0);
t296 = t295-t48.*t291.*(1.0./6.0);
t297 = t78.^2;
t298 = mu.*t297;
t299 = t10-t34-t63+t68+t71-t76;
t300 = t4-t30-t60+t66+t69-t73;
t301 = mu.*t78.*t84;
t302 = lambda.*t28.*t78;
t303 = t302-mu.*t2.*t74;
t304 = t48.*t303.*(1.0./6.0);
t305 = lambda.*t12.*t78;
t306 = t305-mu.*t2.*t72;
t307 = t304-t51.*t52.*t306.*(1.0./6.0);
t308 = lambda.*t32.*t78;
t309 = t308-mu.*t2.*t75;
t310 = lambda.*t6.*t78;
t311 = t310-mu.*t2.*t70;
t312 = t51.*t52.*t311.*(1.0./6.0);
t313 = t312-t48.*t309.*(1.0./6.0);
t314 = t124+t125+t126+t127;
t315 = t51.*t52.*(t128+t129-t2.*t8.*t78).*(1.0./6.0);
t316 = t315-t48.*t314.*(1.0./6.0);
t317 = lambda.*t64.*t78;
t318 = t317-mu.*t51.*t72;
t319 = t51.*t52.*t318.*(1.0./6.0);
t320 = mu.*t48.*t51.*t74.*(1.0./6.0);
t321 = t319+t320;
t322 = lambda.*t61.*t78;
t323 = t322-mu.*t51.*t70;
t324 = t51.*t52.*t323.*(-1.0./6.0)-mu.*t48.*t51.*t75.*(1.0./6.0);
t325 = t204+t205;
t326 = t48.*t325.*(1.0./6.0);
t327 = t206+t207-t8.*t51.*t78;
t328 = t326-t51.*t52.*t327.*(1.0./6.0);
t329 = lambda.*t74.*t78;
t330 = mu.*t74.*t78;
t331 = t329+t330;
t332 = lambda.*t72.*t78;
t333 = mu.*t72.*t78;
t334 = t51.*t52.*(t332+t333).*(1.0./6.0);
t335 = t334-t48.*t331.*(1.0./6.0);
t336 = lambda.*t75.*t78;
t337 = mu.*t75.*t78;
t338 = t336+t337;
t339 = t48.*t338.*(1.0./6.0);
t340 = lambda.*t70.*t78;
t341 = mu.*t70.*t78;
t342 = t340+t341;
t343 = t339-t51.*t52.*t342.*(1.0./6.0);
t344 = mu.*t70.*t75.*2.0;
t345 = mu.*t72.*t74.*2.0;
t346 = t4-t30-t60+t66+t69-t73;
t347 = t10-t34-t63+t68+t71-t76;
t348 = mu.*t70.*t83;
t349 = mu.*t72.*t82;
t350 = mu.*t70.*t80;
t351 = mu.*t72.*t81;
t352 = t8.*t12.*t82;
t490 = mu.*t32.*t80;
t353 = t48.*(t130+t352-t490-t8.*t28.*t81).*(1.0./6.0);
t354 = t8.*t12.*t81;
t420 = mu.*t2.*t84;
t355 = t132+t354-t420;
t356 = t51.*t52.*t355.*(1.0./6.0);
t357 = t353+t356;
t358 = lambda.*t6.*t82;
t359 = mu.*t12.*t83;
t360 = t358+t359-lambda.*t32.*t81-mu.*t28.*t80;
t361 = lambda.*t6.*t81;
t362 = mu.*t12.*t80;
t363 = t361+t362;
t364 = t48.*t360.*(-1.0./6.0)-t51.*t52.*t363.*(1.0./6.0);
t365 = lambda.*t2.*t82;
t366 = mu.*t28.*(t49-t53+t54-t57-t77+t79);
t367 = t365+t366;
t368 = t48.*t367.*(1.0./6.0);
t369 = lambda.*t2.*t81;
t370 = t369-mu.*t12.*t84;
t371 = t51.*t52.*t370.*(1.0./6.0);
t372 = t368+t371;
t373 = t8.*t64.*t82;
t374 = t208+t373;
t375 = t8.*t64.*t81;
t439 = mu.*t51.*t84;
t376 = t210+t375-t439;
t377 = t48.*t374.*(-1.0./6.0)-t51.*t52.*t376.*(1.0./6.0);
t378 = lambda.*t61.*t82;
t379 = mu.*t64.*t83;
t380 = t378+t379;
t381 = t48.*t380.*(1.0./6.0);
t382 = lambda.*t61.*t81;
t383 = mu.*t64.*t80;
t384 = t382+t383;
t385 = t51.*t52.*t384.*(1.0./6.0);
t386 = t381+t385;
t387 = lambda.*t51.*t81;
t388 = t387-mu.*t64.*t84;
t389 = t51.*t52.*t388.*(-1.0./6.0)-lambda.*t48.*t51.*t82.*(1.0./6.0);
t390 = t8.*t72.*t82;
t517 = mu.*t75.*t80;
t391 = t348+t390-t517-t8.*t74.*t81;
t392 = t8.*t72.*t81;
t393 = t301+t350+t392;
t394 = t48.*t391.*(-1.0./6.0)-t51.*t52.*t393.*(1.0./6.0);
t395 = lambda.*t70.*t82;
t396 = mu.*t72.*t83;
t397 = t48.*(t395+t396-lambda.*t75.*t81-mu.*t74.*t80).*(1.0./6.0);
t398 = lambda.*t70.*t81;
t399 = mu.*t72.*t80;
t400 = t51.*t52.*(t398+t399).*(1.0./6.0);
t401 = t397+t400;
t402 = lambda.*t78.*t82;
t403 = t402-mu.*t74.*t84;
t404 = t48.*t403.*(1.0./6.0);
t405 = lambda.*t78.*t81;
t406 = mu.*(t10-t34-t63+t68+t71-t76).*(t49-t53+t54-t57-t77+t79);
t407 = t405+t406;
t408 = t51.*t52.*t407.*(1.0./6.0);
t409 = t404+t408;
t410 = t49-t53+t54-t57-t77+t79;
t411 = mu.*t6.*t82;
t412 = lambda.*t12.*t83;
t413 = t411+t412-lambda.*t28.*t80-mu.*t32.*t81;
t414 = lambda.*t12.*t80;
t415 = mu.*t6.*t81;
t416 = t414+t415;
t417 = t48.*t413.*(-1.0./6.0)-t51.*t52.*t416.*(1.0./6.0);
t418 = t6.*t8.*t83;
t491 = mu.*t28.*t81;
t419 = t48.*(t131+t418-t491-t8.*t32.*t80).*(1.0./6.0);
t421 = t6.*t8.*t80;
t422 = lambda.*t2.*t83;
t423 = mu.*t32.*t84;
t424 = t422+t423;
t425 = lambda.*t2.*t80;
t426 = t425-mu.*t6.*t84;
t427 = t48.*t424.*(-1.0./6.0)-t51.*t52.*t426.*(1.0./6.0);
t428 = mu.*t61.*t82;
t429 = lambda.*t64.*t83;
t430 = t428+t429;
t431 = t48.*t430.*(1.0./6.0);
t432 = lambda.*t64.*t80;
t433 = mu.*t61.*t81;
t434 = t432+t433;
t435 = t51.*t52.*t434.*(1.0./6.0);
t436 = t431+t435;
t437 = t8.*t61.*t83;
t438 = t209+t437;
t440 = t8.*t61.*t80;
t441 = lambda.*t51.*t80;
t442 = t441-mu.*t61.*t84;
t443 = t51.*t52.*t442.*(1.0./6.0);
t444 = lambda.*t48.*t51.*t83.*(1.0./6.0);
t445 = t443+t444;
t446 = lambda.*t74.*t80;
t447 = mu.*t75.*t81;
t448 = t446+t447-lambda.*t72.*t83-mu.*t70.*t82;
t449 = lambda.*t72.*t80;
t450 = mu.*t70.*t81;
t451 = t51.*t52.*(t449+t450).*(1.0./6.0);
t452 = t451-t48.*t448.*(1.0./6.0);
t453 = t8.*t70.*t83;
t518 = mu.*t74.*t81;
t454 = t349+t453-t518-t8.*t75.*t80;
t455 = t8.*t70.*t80;
t456 = t301+t351+t455;
t457 = t48.*t454.*(-1.0./6.0)-t51.*t52.*t456.*(1.0./6.0);
t458 = lambda.*t78.*t83;
t459 = t458-mu.*t75.*t84;
t460 = lambda.*t78.*t80;
t461 = mu.*t70.*t84;
t462 = t460+t461;
t463 = t48.*t459.*(-1.0./6.0)-t51.*t52.*t462.*(1.0./6.0);
t464 = lambda.*t80.*t82;
t465 = mu.*t80.*t82;
t466 = lambda.*t81.*t83;
t467 = mu.*t81.*t83;
t468 = t464+t465+t466+t467;
t469 = t48.*t468.*(1.0./6.0);
t470 = lambda.*t80.*t81;
t471 = mu.*t80.*t81;
t472 = t470+t471;
t473 = t51.*t52.*t472.*(1.0./6.0);
t474 = t469+t473;
t475 = t49-t53+t54-t57-t77+t79;
t476 = t81.^2;
t477 = t80.^2;
t478 = lambda.*t28.*t84;
t479 = mu.*t2.*t82;
t480 = t48.*(t478+t479).*(1.0./6.0);
t481 = lambda.*t12.*t84;
t482 = t481-mu.*t2.*t81;
t483 = t480-t51.*t52.*t482.*(1.0./6.0);
t484 = lambda.*t32.*t84;
t485 = mu.*t2.*t83;
t486 = t484+t485;
t487 = lambda.*t6.*t84;
t488 = t51.*t52.*(t487-mu.*t2.*t80).*(1.0./6.0);
t489 = t488-t48.*t486.*(1.0./6.0);
t492 = t132+t133-t2.*t8.*t84;
t493 = t51.*t52.*t492.*(1.0./6.0);
t494 = lambda.*t64.*t84;
t495 = t51.*t52.*(t494-mu.*t51.*t81).*(1.0./6.0);
t496 = t495-mu.*t48.*t51.*t82.*(1.0./6.0);
t497 = lambda.*t61.*t84;
t498 = t497-mu.*t51.*t80;
t499 = mu.*t48.*t51.*t83.*(1.0./6.0);
t500 = t499-t51.*t52.*t498.*(1.0./6.0);
t501 = t208+t209;
t502 = t210+t211-t8.*t51.*t84;
t503 = t48.*t501.*(-1.0./6.0)-t51.*t52.*t502.*(1.0./6.0);
t504 = lambda.*t74.*t84;
t505 = t504-mu.*t78.*t82;
t506 = lambda.*(t10-t34-t63+t68+t71-t76).*(t49-t53+t54-t57-t77+t79);
t507 = mu.*t78.*t81;
t508 = t506+t507;
t509 = t51.*t52.*t508.*(1.0./6.0);
t510 = t509-t48.*t505.*(1.0./6.0);
t511 = lambda.*t75.*t84;
t512 = t48.*(t511-mu.*t78.*t83).*(1.0./6.0);
t513 = mu.*t78.*t80;
t514 = lambda.*t70.*t84;
t515 = t513+t514;
t516 = t512-t51.*t52.*t515.*(1.0./6.0);
t519 = t8.*t78.*t84;
t520 = t350+t351+t519;
t521 = lambda.*t82.*t84;
t522 = mu.*t82.*t84;
t523 = t48.*(t521+t522).*(1.0./6.0);
t524 = lambda.*t81.*t84;
t525 = mu.*t81.*t84;
t526 = t51.*t52.*(t524+t525).*(1.0./6.0);
t527 = t523+t526;
t528 = lambda.*t83.*t84;
t529 = mu.*t83.*t84;
t530 = t528+t529;
t531 = lambda.*t80.*t84;
t532 = mu.*t80.*t84;
t533 = t531+t532;
t534 = t48.*t530.*(-1.0./6.0)-t51.*t52.*t533.*(1.0./6.0);
t535 = mu.*t80.*t83.*2.0;
t536 = mu.*t81.*t82.*2.0;
t537 = mu.*t477;
t538 = mu.*t476;
t539 = t49-t53+t54-t57-t77+t79;
der_st_dz2 = reshape([t48.*(t116+t8.*t12.*t28.*2.0).*(1.0./6.0)-t51.*t52.*(t96+t118+t8.*t97).*(1.0./6.0),t94,t107,t139,t160,t184,t217,t260,t307,t357,t417,t483,t94,t48.*(t117+t6.*t8.*t32.*2.0).*(1.0./6.0)-t51.*t52.*(t96+t119+t8.*t98).*(1.0./6.0),t115,t147,t166,t189,t227,t48.*t263.*(-1.0./6.0)+t51.*t52.*(t129-t264+t265).*(1.0./6.0),t313,t364,t419+t51.*t52.*(t133-t420+t421).*(1.0./6.0),t489,t107,t115,t48.*(t116+t117).*(1.0./6.0)-t51.*t52.*(t118+t119+t8.*t95).*(1.0./6.0),t152,t171,t194,t232,t271,t316,t372,t427,t493+t48.*(t130+t131-t490-t491).*(1.0./6.0),t139,t147,t152,t51.*t52.*(t177+t202+t8.*t178).*(-1.0./6.0),t175,t51.*t52.*t197.*(-1.0./6.0),t238,t278,t321,t377,t436,t496,t160,t166,t171,t175,t51.*t52.*(t177+t203+t8.*t179).*(-1.0./6.0),t201,t245,t281-t51.*t52.*(t207-t282+t283).*(1.0./6.0),t324,t386,t48.*t438.*(-1.0./6.0)-t51.*t52.*(t211-t439+t440).*(1.0./6.0),t500,t184,t189,t194,t51.*t52.*t197.*(-1.0./6.0),t201,t51.*t52.*(t202+t203+t8.*t176).*(-1.0./6.0),t249,t286,t328,t389,t445,t503,t217,t227,t232,t238,t245,t249,t48.*(t344+t8.*t72.*t74.*2.0).*(1.0./6.0)-t51.*t52.*(t298+mu.*t250.^2+t8.*t251.^2).*(1.0./6.0),t296,t335,t394,t452,t510,t260,t48.*t263.*(-1.0./6.0)+t51.*t52.*(t129+t265-mu.*t2.*t78).*(1.0./6.0),t271,t278,t281-t51.*t52.*(t207+t283-mu.*t51.*t78).*(1.0./6.0),t286,t296,t48.*(t345+t8.*t70.*t75.*2.0).*(1.0./6.0)-t51.*t52.*(t298+mu.*t299.^2+t8.*t300.^2).*(1.0./6.0),t343,t401,t457,t516,t307,t313,t316,t321,t324,t328,t335,t343,t48.*(t344+t345).*(1.0./6.0)-t51.*t52.*(t8.*t297+mu.*t346.^2+mu.*t347.^2).*(1.0./6.0),t409,t463,t48.*(t348+t349-t517-t518).*(-1.0./6.0)-t51.*t52.*t520.*(1.0./6.0),t357,t364,t372,t377,t386,t389,t394,t401,t409,t48.*(t535+t8.*t81.*t82.*2.0).*(-1.0./6.0)-t51.*t52.*(t537+t8.*t476+mu.*t410.^2).*(1.0./6.0),t474,t527,t417,t419+t51.*t52.*(t133+t421-mu.*t2.*t84).*(1.0./6.0),t427,t436,t48.*t438.*(-1.0./6.0)-t51.*t52.*(t211+t440-mu.*t51.*t84).*(1.0./6.0),t445,t452,t457,t463,t474,t48.*(t536+t8.*t80.*t83.*2.0).*(-1.0./6.0)-t51.*t52.*(t538+t8.*t477+mu.*t475.^2).*(1.0./6.0),t534,t483,t489,t493+t48.*(t130+t131-mu.*t28.*t81-mu.*t32.*t80).*(1.0./6.0),t496,t500,t503,t510,t516,t48.*(t348+t349-mu.*t74.*t81-mu.*t75.*t80).*(-1.0./6.0)-t51.*t52.*t520.*(1.0./6.0),t527,t534,t48.*(t535+t536).*(-1.0./6.0)-t51.*t52.*(t537+t538+t8.*t539.^2).*(1.0./6.0)],[12, 12]);
