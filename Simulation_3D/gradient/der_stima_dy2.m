function der_st_dy2 = der_stima_dy2(lambda,mu,V)
%DER_STIMA_DY2
%    DER_ST_DY2 = DER_STIMA_DY2(LAMBDA,MU,X1,X2,X3,X4,Y1,Y2,Y3,Y4,Z1,Z2,Z3,Z4)

%    This function was generated by the Symbolic Math Toolbox version 5.8.
%    02-Jun-2015 16:28:22
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

t2 = x2.*y3;
t3 = x4.*y2;
t4 = x3.*y4;
t53 = x3.*y2;
t54 = x2.*y4;
t55 = x4.*y3;
t5 = t2+t3+t4-t53-t54-t55;
t13 = x3.*z4;
t14 = x4.*z3;
t29 = x2.*z3;
t30 = x3.*z2;
t31 = x2.*z4;
t32 = x4.*z2;
t6 = t13-t14+t29-t30-t31+t32;
t7 = mu.*2.0;
t8 = lambda+t7;
t9 = y2.*z3;
t10 = y4.*z2;
t11 = y3.*z4;
t46 = y3.*z2;
t47 = y2.*z4;
t48 = y4.*z3;
t12 = t9+t10+t11-t46-t47-t48;
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
t33 = x1.*y3.*z2;
t34 = x2.*y1.*z3;
t35 = x3.*y2.*z1;
t36 = x1.*y2.*z4;
t37 = x2.*y4.*z1;
t38 = x4.*y1.*z2;
t39 = x1.*y4.*z3;
t40 = x3.*y1.*z4;
t41 = x4.*y3.*z1;
t42 = x2.*y3.*z4;
t43 = x3.*y4.*z2;
t44 = x4.*y2.*z3;
t27 = t15+t16+t17+t18+t19+t20+t21+t22+t23+t24+t25+t26-t33-t34-t35-t36-t37-t38-t39-t40-t41-t42-t43-t44;
t28 = z3-z4;
t45 = 1.0./t27;
t49 = x1.*z3;
t50 = x4.*z1;
t57 = x3.*z1;
t58 = x1.*z4;
t51 = t13-t14+t49+t50-t57-t58;
t52 = 1.0./t27.^2;
t56 = x3-x4;
t59 = x1.*y3;
t60 = x4.*y1;
t67 = x3.*y1;
t68 = x1.*y4;
t61 = t4-t55+t59+t60-t67-t68;
t62 = y1.*z3;
t63 = y4.*z1;
t65 = y3.*z1;
t66 = y1.*z4;
t64 = t11-t48+t62+t63-t65-t66;
t69 = x2.*y1;
t78 = x1.*y2;
t70 = t3-t54-t60+t68+t69-t78;
t71 = y2.*z1;
t77 = y1.*z2;
t72 = t10-t47-t63+t66+t71-t77;
t73 = x1.*z2;
t76 = x2.*z1;
t74 = t31-t32+t50-t58+t73-t76;
t75 = z1-z4;
t79 = x1-x4;
t80 = t2-t53-t59+t67-t69+t78;
t81 = t9-t46-t62+t65-t71+t77;
t82 = t29-t30-t49+t57+t73-t76;
t83 = z1-z3;
t84 = x1-x3;
t85 = lambda.*t6.*t28;
t86 = mu.*t6.*t28;
t87 = t85+t86;
t88 = t45.*t87.*(1.0./6.0);
t89 = lambda.*t6.*t12;
t90 = mu.*t6.*t12;
t91 = t89+t90;
t92 = t88-t51.*t52.*t91.*(1.0./6.0);
t93 = mu.*t5.*t56.*2.0;
t94 = t5.^2;
t95 = mu.*t94;
t96 = t12.^2;
t97 = t6.^2;
t98 = mu.*t56.*t61;
t99 = mu.*t5.*t61;
t100 = mu.*t56.*t70;
t101 = mu.*t5.*t79;
t102 = mu.*t5.*t70;
t103 = mu.*t5.*t84;
t104 = mu.*t5.*t80;
t105 = lambda.*t5.*t28;
t106 = mu.*t5.*t28;
t107 = t105+t106-lambda.*t12.*t56-mu.*t12.*t56;
t108 = lambda.*t5.*t12;
t109 = mu.*t5.*t12;
t110 = t108+t109;
t111 = t51.*t52.*t110.*(1.0./6.0);
t112 = t111-t45.*t107.*(1.0./6.0);
t113 = lambda.*t6.*t56;
t114 = mu.*t6.*t56;
t115 = t113+t114;
t116 = lambda.*t5.*t6;
t117 = mu.*t5.*t6;
t118 = t116+t117;
t119 = t45.*t115.*(-1.0./6.0)-t51.*t52.*t118.*(1.0./6.0);
t120 = mu.*t97;
t121 = mu.*t96;
t122 = mu.*t6.*t51;
t123 = mu.*t12.*t64;
t124 = mu.*t12.*t72;
t125 = mu.*t28.*t81;
t126 = mu.*t6.*t82;
t127 = mu.*t12.*t81;
t128 = t98-t8.*t28.*t64;
t129 = t8.*t12.*t64;
t130 = t99+t122+t129;
t131 = t45.*t128.*(-1.0./6.0)-t51.*t52.*t130.*(1.0./6.0);
t132 = lambda.*t6.*t64;
t133 = mu.*t12.*t51;
t134 = t132+t133;
t135 = t51.*t52.*t134.*(1.0./6.0);
t136 = t135-mu.*t28.*t45.*t51.*(1.0./6.0);
t137 = lambda.*t56.*t64;
t138 = t137-mu.*t28.*t61;
t139 = lambda.*t5.*t64;
t140 = mu.*t12.*t61;
t141 = t139+t140;
t142 = t45.*t138.*(-1.0./6.0)-t51.*t52.*t141.*(1.0./6.0);
t143 = lambda.*t12.*t51;
t144 = mu.*t6.*t64;
t145 = t143+t144;
t146 = t51.*t52.*t145.*(1.0./6.0);
t147 = t146-lambda.*t28.*t45.*t51.*(1.0./6.0);
t148 = mu.*t28.*t64;
t149 = t6.*t8.*t51;
t150 = t99+t123+t149;
t151 = lambda.*t5.*t51;
t152 = mu.*t6.*t61;
t153 = t151+t152;
t154 = t51.*t52.*t153.*(1.0./6.0);
t155 = lambda.*t45.*t51.*t56.*(1.0./6.0);
t156 = t154+t155;
t157 = lambda.*t51.*t64;
t158 = mu.*t51.*t64;
t159 = t157+t158;
t160 = t61.^2;
t161 = mu.*t160;
t162 = t64.^2;
t163 = t51.^2;
t164 = mu.*t61.*t79;
t165 = mu.*t61.*t70;
t166 = mu.*t61.*t84;
t167 = mu.*t61.*t80;
t168 = lambda.*t28.*t61;
t169 = t168-mu.*t56.*t64;
t170 = t45.*t169.*(1.0./6.0);
t171 = lambda.*t12.*t61;
t172 = mu.*t5.*t64;
t173 = t171+t172;
t174 = t170-t51.*t52.*t173.*(1.0./6.0);
t175 = lambda.*t6.*t61;
t176 = mu.*t5.*t51;
t177 = t175+t176;
t178 = t51.*t52.*t177.*(1.0./6.0);
t179 = mu.*t45.*t51.*t56.*(1.0./6.0);
t180 = t178+t179;
t181 = t148-t8.*t56.*t61;
t182 = t45.*t181.*(1.0./6.0);
t183 = t5.*t8.*t61;
t184 = t122+t123+t183;
t185 = t182-t51.*t52.*t184.*(1.0./6.0);
t186 = lambda.*t61.*t64;
t187 = mu.*t61.*t64;
t188 = t186+t187;
t189 = t51.*t52.*t188.*(1.0./6.0);
t190 = lambda.*t51.*t61;
t191 = mu.*t51.*t61;
t192 = t190+t191;
t193 = mu.*t163;
t194 = mu.*t162;
t195 = mu.*t64.*t72;
t196 = mu.*t51.*t82;
t197 = mu.*t64.*t81;
t198 = t100+t101-t8.*t12.*t75-t8.*t28.*t72;
t199 = t8.*t12.*t72;
t283 = mu.*t6.*t74;
t200 = t102+t199-t283;
t201 = t45.*t198.*(-1.0./6.0)-t51.*t52.*t200.*(1.0./6.0);
t202 = lambda.*t6.*t75;
t203 = t202-mu.*t28.*t74;
t204 = lambda.*t6.*t72;
t205 = t51.*t52.*(t204-mu.*t12.*t74).*(1.0./6.0);
t206 = t205-t45.*t203.*(1.0./6.0);
t207 = lambda.*t5.*t75;
t208 = mu.*t28.*(t3-t54-t60+t68+t69-t78);
t209 = t207+t208-lambda.*t56.*t72-mu.*t12.*t79;
t210 = t45.*t209.*(1.0./6.0);
t211 = lambda.*t5.*t72;
t212 = mu.*t12.*t70;
t213 = t211+t212;
t214 = t210-t51.*t52.*t213.*(1.0./6.0);
t215 = t164-t8.*t64.*t75;
t216 = t45.*t215.*(1.0./6.0);
t217 = t8.*t64.*t72;
t295 = mu.*t51.*t74;
t218 = t51.*t52.*(t165+t217-t295).*(1.0./6.0);
t219 = t216+t218;
t220 = lambda.*t51.*t72;
t221 = t220-mu.*t64.*t74;
t222 = lambda.*t45.*t51.*t75.*(1.0./6.0);
t223 = t222-t51.*t52.*t221.*(1.0./6.0);
t224 = lambda.*t61.*t75;
t225 = t224-mu.*t64.*t79;
t226 = lambda.*t61.*t72;
t227 = mu.*t64.*t70;
t228 = t51.*t52.*(t226+t227).*(1.0./6.0);
t229 = t228-t45.*t225.*(1.0./6.0);
t230 = t3-t54-t60+t68+t69-t78;
t231 = t10-t47-t63+t66+t71-t77;
t232 = lambda.*t28.*t74;
t233 = t232-mu.*t6.*t75;
t234 = t45.*t233.*(1.0./6.0);
t235 = lambda.*t12.*t74;
t236 = t235-mu.*t6.*t72;
t237 = t234-t51.*t52.*t236.*(1.0./6.0);
t238 = mu.*t12.*t75;
t239 = mu.*t28.*t72;
t240 = t102+t124-t6.*t8.*t74;
t241 = lambda.*t56.*t74;
t242 = t241-mu.*t6.*t79;
t243 = lambda.*t5.*t74;
t244 = t243-mu.*t6.*t70;
t245 = t45.*t242.*(-1.0./6.0)-t51.*t52.*t244.*(1.0./6.0);
t246 = lambda.*t64.*t74;
t247 = t246-mu.*t51.*t72;
t248 = t51.*t52.*t247.*(1.0./6.0);
t249 = mu.*t45.*t51.*t75.*(1.0./6.0);
t250 = t248+t249;
t251 = mu.*t64.*t75;
t252 = t51.*t52.*(t165+t195-t8.*t51.*t74).*(1.0./6.0);
t253 = lambda.*t61.*t74;
t254 = t253-mu.*t51.*t70;
t255 = t51.*t52.*t254.*(1.0./6.0);
t256 = t255-mu.*t45.*t51.*t79.*(1.0./6.0);
t257 = lambda.*t74.*t75;
t258 = mu.*t74.*t75;
t259 = t257+t258;
t260 = lambda.*t72.*t74;
t261 = mu.*t72.*t74;
t262 = t51.*t52.*(t260+t261).*(1.0./6.0);
t263 = t262-t45.*t259.*(1.0./6.0);
t264 = mu.*t70.*t79.*2.0;
t265 = t3-t54-t60+t68+t69-t78;
t266 = t10-t47-t63+t66+t71-t77;
t267 = t74.^2;
t268 = mu.*t79.*t80;
t269 = mu.*t70.*t80;
t270 = lambda.*t12.*t79;
t271 = mu.*t56.*t72;
t272 = t270+t271-lambda.*t28.*t70-mu.*t5.*t75;
t273 = lambda.*t12.*t70;
t274 = mu.*t5.*t72;
t275 = t273+t274;
t276 = t45.*t272.*(-1.0./6.0)-t51.*t52.*t275.*(1.0./6.0);
t277 = lambda.*t6.*t79;
t278 = t277-mu.*t56.*t74;
t279 = t45.*t278.*(1.0./6.0);
t280 = lambda.*t6.*t70;
t281 = t51.*t52.*(t280-mu.*t5.*t74).*(1.0./6.0);
t282 = t279+t281;
t284 = t5.*t8.*t70;
t285 = lambda.*t64.*t79;
t286 = t45.*(t285-mu.*t61.*t75).*(1.0./6.0);
t287 = lambda.*t64.*t70;
t288 = mu.*t61.*t72;
t289 = t51.*t52.*(t287+t288).*(1.0./6.0);
t290 = t286+t289;
t291 = lambda.*t51.*t70;
t292 = t291-mu.*t61.*t74;
t293 = t51.*t52.*t292.*(-1.0./6.0)-lambda.*t45.*t51.*t79.*(1.0./6.0);
t294 = t251-t8.*t61.*t79;
t296 = t8.*t61.*t70;
t297 = lambda.*t79.*(t10-t47-t63+t66+t71-t77);
t298 = mu.*t79.*(t10-t47-t63+t66+t71-t77);
t299 = t297+t298-lambda.*t70.*t75-mu.*t70.*t75;
t300 = t45.*t299.*(1.0./6.0);
t301 = lambda.*(t3-t54-t60+t68+t69-t78).*(t10-t47-t63+t66+t71-t77);
t302 = mu.*(t3-t54-t60+t68+t69-t78).*(t10-t47-t63+t66+t71-t77);
t303 = t301+t302;
t304 = t51.*t52.*t303.*(1.0./6.0);
t305 = t300+t304;
t306 = lambda.*t74.*t79;
t307 = mu.*t74.*t79;
t308 = t306+t307;
t309 = t45.*t308.*(1.0./6.0);
t310 = lambda.*t70.*t74;
t311 = mu.*t70.*t74;
t312 = t51.*t52.*(t310+t311).*(1.0./6.0);
t313 = t309+t312;
t314 = mu.*t267;
t315 = t10-t47-t63+t66+t71-t77;
t316 = t3-t54-t60+t68+t69-t78;
t317 = mu.*t72.*t81;
t318 = t8.*t28.*t81;
t381 = mu.*t56.*t80;
t319 = t45.*(t103+t318-t381-t8.*t12.*t83).*(1.0./6.0);
t320 = t8.*t12.*t81;
t321 = t104+t126+t320;
t322 = t319-t51.*t52.*t321.*(1.0./6.0);
t323 = lambda.*t6.*t83;
t324 = t323-mu.*t28.*t82;
t325 = t45.*t324.*(1.0./6.0);
t326 = lambda.*t6.*t81;
t327 = mu.*t12.*t82;
t328 = t326+t327;
t329 = t51.*t52.*t328.*(1.0./6.0);
t330 = t325+t329;
t331 = lambda.*t5.*t83;
t332 = lambda.*t56.*t81;
t333 = t331+t332-mu.*t12.*t84-mu.*t28.*t80;
t334 = lambda.*t5.*t81;
t335 = mu.*t12.*t80;
t336 = t334+t335;
t337 = t45.*t333.*(-1.0./6.0)-t51.*t52.*t336.*(1.0./6.0);
t338 = t166-t8.*t64.*t83;
t339 = t8.*t64.*t81;
t340 = t167+t196+t339;
t341 = t51.*t52.*t340.*(1.0./6.0);
t342 = t341-t45.*t338.*(1.0./6.0);
t343 = lambda.*t51.*t81;
t344 = mu.*t64.*t82;
t345 = t343+t344;
t346 = t51.*t52.*t345.*(-1.0./6.0)-lambda.*t45.*t51.*t83.*(1.0./6.0);
t347 = lambda.*t61.*t83;
t348 = t347-mu.*t64.*t84;
t349 = t45.*t348.*(1.0./6.0);
t350 = lambda.*t61.*t81;
t351 = mu.*t64.*t80;
t352 = t350+t351;
t353 = t51.*t52.*t352.*(1.0./6.0);
t354 = t349+t353;
t355 = t8.*t83.*(t10-t47-t63+t66+t71-t77);
t412 = mu.*t70.*t84;
t356 = t268+t355-t412-t8.*t75.*t81;
t357 = t45.*t356.*(1.0./6.0);
t358 = t8.*t72.*t81;
t488 = mu.*t74.*t82;
t359 = t51.*t52.*(t269+t358-t488).*(1.0./6.0);
t360 = t357+t359;
t361 = lambda.*t74.*t83;
t362 = mu.*t75.*t82;
t363 = t361+t362;
t364 = t45.*t363.*(1.0./6.0);
t365 = lambda.*t74.*t81;
t366 = t365-mu.*t72.*t82;
t367 = t51.*t52.*t366.*(1.0./6.0);
t368 = t364+t367;
t369 = lambda.*t70.*t83;
t370 = lambda.*t79.*t81;
t371 = lambda.*t70.*t81;
t372 = mu.*t72.*t80;
t373 = t51.*t52.*(t371+t372).*(1.0./6.0);
t374 = lambda.*t28.*t82;
t375 = t374-mu.*t6.*t83;
t376 = lambda.*t12.*t82;
t377 = mu.*t6.*t81;
t378 = t376+t377;
t379 = t51.*t52.*t378.*(1.0./6.0);
t380 = t379-t45.*t375.*(1.0./6.0);
t382 = t6.*t8.*t82;
t383 = t104+t127+t382;
t384 = lambda.*t56.*t82;
t385 = t384-mu.*t6.*t84;
t386 = t45.*t385.*(1.0./6.0);
t387 = lambda.*t5.*t82;
t388 = mu.*t6.*t80;
t389 = t387+t388;
t390 = t51.*t52.*t389.*(1.0./6.0);
t391 = t386+t390;
t392 = lambda.*t64.*t82;
t393 = mu.*t51.*t81;
t394 = t392+t393;
t395 = t51.*t52.*t394.*(-1.0./6.0)-mu.*t45.*t51.*t83.*(1.0./6.0);
t396 = mu.*t64.*t83;
t397 = t8.*t51.*t82;
t398 = t167+t197+t397;
t399 = t51.*t52.*t398.*(1.0./6.0);
t400 = lambda.*t61.*t82;
t401 = mu.*t51.*t80;
t402 = t400+t401;
t403 = mu.*t45.*t51.*t84.*(1.0./6.0);
t404 = t403-t51.*t52.*t402.*(1.0./6.0);
t405 = lambda.*t75.*t82;
t406 = mu.*t74.*t83;
t407 = t405+t406;
t408 = t45.*t407.*(1.0./6.0);
t409 = lambda.*t72.*t82;
t410 = t409-mu.*t74.*t81;
t411 = t408-t51.*t52.*t410.*(1.0./6.0);
t413 = mu.*t83.*(t10-t47-t63+t66+t71-t77);
t414 = t51.*t52.*(t269+t317-t8.*t74.*t82).*(1.0./6.0);
t415 = lambda.*t79.*t82;
t416 = mu.*t74.*t84;
t417 = t415+t416;
t418 = lambda.*t70.*t82;
t419 = t418-mu.*t74.*t80;
t420 = t45.*t417.*(-1.0./6.0)-t51.*t52.*t419.*(1.0./6.0);
t421 = lambda.*t82.*t83;
t422 = mu.*t82.*t83;
t423 = t421+t422;
t424 = lambda.*t81.*t82;
t425 = mu.*t81.*t82;
t426 = t424+t425;
t427 = t45.*t423.*(-1.0./6.0)-t51.*t52.*t426.*(1.0./6.0);
t428 = mu.*t80.*t84.*2.0;
t429 = t80.^2;
t430 = mu.*t429;
t431 = t81.^2;
t432 = t82.^2;
t433 = lambda.*t28.*t80;
t434 = lambda.*t12.*t84;
t435 = t433+t434-mu.*t5.*t83-mu.*t56.*t81;
t436 = t45.*t435.*(1.0./6.0);
t437 = lambda.*t12.*t80;
t438 = mu.*t5.*t81;
t439 = t437+t438;
t440 = t436-t51.*t52.*t439.*(1.0./6.0);
t441 = lambda.*t6.*t84;
t442 = t441-mu.*t56.*t82;
t443 = lambda.*t6.*t80;
t444 = mu.*t5.*t82;
t445 = t443+t444;
t446 = t51.*t52.*t445.*(1.0./6.0);
t447 = t446-t45.*t442.*(1.0./6.0);
t448 = t5.*t8.*t84;
t449 = t125+t448-mu.*t12.*t83-t8.*t56.*t80;
t450 = t45.*t449.*(1.0./6.0);
t451 = t5.*t8.*t80;
t452 = t126+t127+t451;
t453 = t450-t51.*t52.*t452.*(1.0./6.0);
t454 = lambda.*t64.*t84;
t455 = t454-mu.*t61.*t83;
t456 = lambda.*t64.*t80;
t457 = mu.*t61.*t81;
t458 = t456+t457;
t459 = t51.*t52.*t458.*(1.0./6.0);
t460 = t459-t45.*t455.*(1.0./6.0);
t461 = lambda.*t51.*t80;
t462 = mu.*t61.*t82;
t463 = t461+t462;
t464 = lambda.*t45.*t51.*t84.*(1.0./6.0);
t465 = t464-t51.*t52.*t463.*(1.0./6.0);
t466 = t396-t8.*t61.*t84;
t467 = t45.*t466.*(1.0./6.0);
t468 = t8.*t61.*t80;
t469 = t196+t197+t468;
t470 = t51.*t52.*t469.*(1.0./6.0);
t471 = t467+t470;
t472 = lambda.*t75.*t80;
t473 = lambda.*t72.*t84;
t474 = t472+t473-mu.*t70.*t83-mu.*t79.*t81;
t475 = lambda.*t72.*t80;
t476 = mu.*t70.*t81;
t477 = t51.*t52.*(t475+t476).*(1.0./6.0);
t478 = t477-t45.*t474.*(1.0./6.0);
t479 = lambda.*t74.*t84;
t480 = mu.*t79.*t82;
t481 = t479+t480;
t482 = lambda.*t74.*t80;
t483 = t482-mu.*t70.*t82;
t484 = t51.*t52.*t483.*(1.0./6.0);
t485 = t484-t45.*t481.*(1.0./6.0);
t486 = mu.*t72.*t83;
t487 = t8.*t79.*t80;
t489 = t8.*t70.*t80;
t490 = lambda.*t80.*t83;
t491 = mu.*t80.*t83;
t492 = t490+t491-lambda.*t81.*t84-mu.*t81.*t84;
t493 = t45.*t492.*(1.0./6.0);
t494 = lambda.*t80.*t81;
t495 = mu.*t80.*t81;
t496 = t494+t495;
t497 = t51.*t52.*t496.*(1.0./6.0);
t498 = t493+t497;
t499 = lambda.*t82.*t84;
t500 = mu.*t82.*t84;
t501 = t499+t500;
t502 = t45.*t501.*(1.0./6.0);
t503 = lambda.*t80.*t82;
t504 = mu.*t80.*t82;
t505 = t503+t504;
t506 = t502-t51.*t52.*t505.*(1.0./6.0);
t507 = mu.*t432;
t508 = mu.*t431;
der_st_dy2 = reshape([t45.*(t93-t8.*t12.*t28.*2.0).*(1.0./6.0)+t51.*t52.*(t95+t120+t8.*t96).*(1.0./6.0),t92,t112,t131,t147,t174,t201,t237,t276,t322,t380,t440,t92,t45.*(t93-mu.*t12.*t28.*2.0).*(1.0./6.0)+t51.*t52.*(t95+t121+t8.*t97).*(1.0./6.0),t119,t136,t45.*(t98-t148).*(-1.0./6.0)-t51.*t52.*t150.*(1.0./6.0),t180,t206,t45.*(t100+t101-t238-t239).*(-1.0./6.0)-t51.*t52.*t240.*(1.0./6.0),t282,t330,t45.*(t103+t125-t381-mu.*t12.*t83).*(1.0./6.0)-t51.*t52.*t383.*(1.0./6.0),t447,t112,t119,t45.*(mu.*t12.*t28.*2.0-t5.*t8.*t56.*2.0).*(-1.0./6.0)+t51.*t52.*(t120+t121+t8.*t94).*(1.0./6.0),t142,t156,t185,t214,t245,t45.*(t238+t239-t8.*t56.*(t3-t54-t60+t68+t69-t78)-t5.*t8.*t79).*(1.0./6.0)-t51.*t52.*(t124-t283+t284).*(1.0./6.0),t337,t391,t453,t131,t136,t142,t51.*t52.*(t161+t193+t8.*t162).*(1.0./6.0),t51.*t52.*t159.*(-1.0./6.0),t189,t219,t250,t290,t342,t395,t460,t147,t45.*(t98-mu.*t28.*t64).*(-1.0./6.0)-t51.*t52.*t150.*(1.0./6.0),t156,t51.*t52.*t159.*(-1.0./6.0),t51.*t52.*(t161+t194+t8.*t163).*(1.0./6.0),t51.*t52.*t192.*(-1.0./6.0),t223,t252+t45.*(t164-t251).*(1.0./6.0),t293,t346,t399-t45.*(t166-t396).*(1.0./6.0),t465,t174,t180,t185,t189,t51.*t52.*t192.*(-1.0./6.0),t51.*t52.*(t193+t194+t8.*t160).*(1.0./6.0),t229,t256,t45.*t294.*(-1.0./6.0)+t51.*t52.*(t195-t295+t296).*(1.0./6.0),t354,t404,t471,t201,t206,t214,t219,t223,t229,t45.*(t264-t8.*t72.*t75.*2.0).*(1.0./6.0)+t51.*t52.*(t314+mu.*t230.^2+t8.*t231.^2).*(1.0./6.0),t263,t305,t360,t411,t478,t237,t45.*(t100+t101-mu.*t12.*t75-mu.*t28.*t72).*(-1.0./6.0)-t51.*t52.*t240.*(1.0./6.0),t245,t250,t252+t45.*(t164-mu.*t64.*t75).*(1.0./6.0),t256,t263,t45.*(t264-mu.*t75.*(t10-t47-t63+t66+t71-t77).*2.0).*(1.0./6.0)+t51.*t52.*(t8.*t267+mu.*t265.^2+mu.*t266.^2).*(1.0./6.0),t313,t368,t414+t45.*(t268-t412+t413-mu.*t75.*t81).*(1.0./6.0),t485,t276,t282,t45.*(t238+t239-t5.*t8.*t79-t8.*t56.*t70).*(1.0./6.0)-t51.*t52.*(t124+t284-mu.*t6.*t74).*(1.0./6.0),t290,t293,t45.*t294.*(-1.0./6.0)+t51.*t52.*(t195+t296-mu.*t51.*t74).*(1.0./6.0),t305,t313,t45.*(t8.*t79.*(t3-t54-t60+t68+t69-t78).*2.0-mu.*t72.*t75.*2.0).*(1.0./6.0)+t51.*t52.*(t314+mu.*t315.^2+t8.*t316.^2).*(1.0./6.0),t373+t45.*(t369+t370-mu.*t84.*(t10-t47-t63+t66+t71-t77)-mu.*t75.*t80).*(1.0./6.0),t420,t45.*(t486+t487-t8.*t84.*(t3-t54-t60+t68+t69-t78)-mu.*t75.*t81).*(1.0./6.0)+t51.*t52.*(t317-t488+t489).*(1.0./6.0),t322,t330,t337,t342,t346,t354,t360,t368,t373+t45.*(t369+t370-mu.*t75.*t80-mu.*t72.*t84).*(1.0./6.0),t45.*(t428-t8.*t81.*t83.*2.0).*(-1.0./6.0)+t51.*t52.*(t430+t507+t8.*t431).*(1.0./6.0),t427,t498,t380,t45.*(t103+t125-mu.*t12.*t83-mu.*t56.*t80).*(1.0./6.0)-t51.*t52.*t383.*(1.0./6.0),t391,t395,t399-t45.*(t166-mu.*t64.*t83).*(1.0./6.0),t404,t411,t414+t45.*(t268+t413-mu.*t70.*t84-mu.*t75.*t81).*(1.0./6.0),t420,t427,t45.*(t428-mu.*t81.*t83.*2.0).*(-1.0./6.0)+t51.*t52.*(t430+t508+t8.*t432).*(1.0./6.0),t506,t440,t447,t453,t460,t465,t471,t478,t485,t45.*(t486+t487-mu.*t75.*t81-t8.*t70.*t84).*(1.0./6.0)+t51.*t52.*(t317+t489-mu.*t74.*t82).*(1.0./6.0),t498,t506,t45.*(mu.*t81.*t83.*2.0-t8.*t80.*t84.*2.0).*(1.0./6.0)+t51.*t52.*(t507+t508+t8.*t429).*(1.0./6.0)],[12, 12]);
